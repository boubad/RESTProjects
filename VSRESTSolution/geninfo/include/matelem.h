#pragma once
#ifndef __MATELEM_H__
#define __MATELEM_H__
/////////////////////////////
#include "matdata.h"
#include "mattree.h"
#include "matelem_result.h"
//////////////////////
namespace info {
	//////////////////////////////////
	template <typename INDEXTYPE, typename DISTANCETYPE, typename CRITERIATYPE>
	class MatElem : CancellableObject {
	public:
		using index_type = INDEXTYPE;
		using criteria_type = CRITERIATYPE;
		using distance_type = DISTANCETYPE;
		using sizets_vector = std::vector<index_type>;
		using distancemap_type = DistanceMap<index_type, distance_type>;
		using cancellableflag_type = CancellableObject::cancellableflag_type;
		//
		using pair_type = std::pair<size_t, size_t>;
		using pairs_type_vector = std::queue<pair_type>;
		using matelemresult_type = MatElemResult<index_type, criteria_type>;
		using matelemresult_type_ptr = std::shared_ptr<matelemresult_type>;
		using task_type = std::future<matelemresult_type_ptr>;
		using matelem_type = MatElem<index_type, distance_type, criteria_type>;
	private:
		distancemap_type *m_pdist;
		matelemresult_type m_res;
		bool m_bnotify;
	protected:
		MatElem(const matelem_type &other) : CancellableObject(other),m_pdist(other.m_pdist), m_res(other.m_res), m_bnotify(other.m_bnotify) {}
		matelem_type & operator=(const matelem_type &other) {
			if (this != &other) {
				CancellableObject::operator=(other);
				this->m_pdist = other.m_pdist;
				this->m_res = other.m_res;
				this->m_bnotify = other.m_bnotify;
			}
		}
		bool find_best_permutations(pairs_type_vector &v, criteria_type &bestCrit) {
			sizets_vector inds;
			this->m_res.indexes(inds);
			const size_t n = inds.size();
			const distancemap_type *pMap = this->m_pdist;
			for (size_t i = 0; i < n; ++i) {
				if (this->is_cancellation_requested()) {
					return (false);
				}
				for (size_t j = 0; j < i; ++j) {
					if (this->is_cancellation_requested()) {
						return (false);
					}
					sizets_vector xinds = inds;
					index_type t = xinds[i];
					xinds[i] = xinds[j];
					xinds[j] = t;
					criteria_type c = 0;
					pMap->compute_criteria(xinds, c);
					if (c < bestCrit) {
						while (!v.empty()) {
							v.pop();
						}
						v.push(std::make_pair(j, i));
						bestCrit = c;
					}
					else if ((c == bestCrit) && (!v.empty())) {
						v.push(std::make_pair(j, i));
					}
				}// j
			}// i
			if (this->is_cancellation_requested()) {
				return (false);
			}
			return (!v.empty());
		}//find_best_permutations
	public:
		template <typename FLOATTYPE, typename STRINGTYPE>
		MatElem(MatData<index_type, FLOATTYPE, distance_type, STRINGTYPE> *pMat,
			DispositionType disp = DispositionType::dispRow,
			cancellableflag_type *pf = nullptr,
			Backgrounder *pBack = nullptr) :CancellableObject(pf, pBack),m_pdist(nullptr),
			m_bnotify(false) {
			assert(pMat != nullptr);
			assert(pMat->is_valid());
			if (disp == DispositionType::dispRow) {
				this->m_pdist = pMat->get_rows_distances_map();
			}
			else {
				this->m_pdist = pMat->get_cols_distances_map();
			}
			sizets_vector oinds;
			this->m_pdist->get_indexes(oinds);
			criteria_type c = 0;
			this->m_pdist->compute_criteria(oinds, c);
			this->m_res.set(c, oinds);
		}// MatElem
		MatElem(distancemap_type *pMap, DispositionType disp = DispositionType::dispUnknown,
			cancellableflag_type *pf = nullptr, Backgrounder *pBack = nullptr) :
			CancellableObject(pf, pBack),m_pdist(pMap),
			m_bnotify(false) {
			assert(pMap != nullptr);
			assert(pMap->is_valid());
			sizets_vector oinds;
			pMap->get_indexes(oinds);
			criteria_type c = 0;
			pMap->compute_criteria(oinds, c);
			this->m_res.set(c, oinds);
		}// MatElem
		virtual ~MatElem() {}
	public:
		DispositionType disposition(void) const {
			return (this->m_res.disposition());
		}
		criteria_type criteria(void) const {
			return (this->m_res.criteria());
		}
		template <typename XU>
		void indexes(std::vector<XU> &oInds) const {
			this->m_res.indexes(oInds);
		}// indexes
		template <typename U, typename X>
		void set(const X &crit, const std::vector<U> &oInds) {
			this->m_res.set(crit, oInds);
		}// set
		bool is_valid(void) const {
			return (this->m_pdist != nullptr);
		}
		matelemresult_type_ptr get_result(void) {
			return (std::make_shared<matelemresult_type>(this->m_res));
		}// get_result
		matelemresult_type_ptr arrange(void) {
			bool done = false;
			matelemresult_type oRet;
			distancemap_type *pMap = this->m_pdist;
			DispositionType disp = this->disposition();
			cancellableflag_type *pf = this->get_cancellable_flag();
			Backgrounder *pb = this->get_backgrounder();
			criteria_type bestCrit = this->criteria();
			while (!done) {
				if (this->is_cancellation_requested()) {
					done = true;
					break;
				}
				pairs_type_vector v;
				if (!this->find_best_permutations(v, bestCrit))
				{
					done = true;
					break;
				}
				pair_type p = v.front();
				v.pop();
				this->m_res.permute(p.first, p.second);
				this->m_res.criteria(bestCrit);
				if (!v.empty()) {
					sizets_vector inds;
					this->m_res.indexes(inds);
					matelemresult_type rbest;
					while (!v.empty()) {
						pair_type pp = v.front();
						v.pop();
						matelem_type xMat(pMap, disp, pf, pb);
						xMat.m_res.set(bestCrit, inds);
						xMat.m_res.permute(pp.first, pp.second);
						matelemresult_type_ptr r = xMat.arrange();
						matelemresult_type *px = r.get();
						if (px != nullptr) {
							rbest += *px;
						}
					}// not empty
					criteria_type cx = rbest.criteria();
					if ((cx > 0) && (cx < this->criteria())) {
						this->m_res = rbest;
					}
				}// not empty
			}// not done
			if (this->is_cancellation_requested()) {
				return (matelemresult_type_ptr());
			}
			return (this->get_result());
		}// arrange

		static matelemresult_type_ptr st_arrange(const distancemap_type *pMap,
			DispositionType disp = DispositionType::dispUnknown,
			cancellableflag_type *pf = nullptr, Backgrounder *pBack = nullptr) {
			matelem_type oMat(pMap, disp, pf, pBack);
			return oMat.arrange();
		}// create
		template <typename FLOATTYPE, typename STRINGTYPE>
		static matelemresult_type_ptr st_arrange_matdata(const MatData<index_type, FLOATTYPE, distance_type, STRINGTYPE> *pMat,
			DispositionType disp = DispositionType::dispRow,
			cancellableflag_type *pf = nullptr, Backgrounder *pBack = nullptr) {
			matelem_type oMat(pMat, disp, pf, pBack);
			return oMat.arrange();
		}// create
		template <typename FLOATTYPE, typename STRINGTYPE>
		static matelemresult_type_ptr st_arrange_hierar(MatData<index_type, FLOATTYPE, distance_type, STRINGTYPE> *pMat,
			DispositionType disp = DispositionType::dispRow,
			cancellableflag_type *pf = nullptr, Backgrounder *pBack = nullptr) {
			MatTree<index_type, FLOATTYPE, distance_type, STRINGTYPE> oTree(pMat, disp, pf,pBack);
			if (!oTree.aggreg()) {
				return (matelemresult_type_ptr());
			}
			sizets_vector oinds;
			oTree.get_ids(oinds);
			matelem_type oMat(pMat, disp, pf, pBack);
			criteria_type c = oMat.m_res.criteria();
			oMat.m_res.set(c, oinds);
			return oMat.arrange();
		}// create
	}; // class MatElem
	   ///////////////////////////////////
}// namespace info
 ///////////////////////////////
#endif // __MATELEM_H__
