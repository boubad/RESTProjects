#pragma once
#ifndef __MATORD_H__
#define __MATORD_H__
/////////////////////////////
#include "matelem.h"
//////////////////////
namespace info {
	//////////////////////////////////
	template <typename INDEXTYPE, typename DISTANCETYPE, typename CRITERIATYPE>
	class MatOrd : public CancellableObject {
	public:
		using index_type = INDEXTYPE;
		using criteria_type = CRITERIATYPE;
		using distance_type = DISTANCETYPE;
		using cancellableflag_type = CancellableObject::cancellableflag_type;
		//
		using sizets_vector = std::vector<index_type>;
		using matelemresult_type = MatElemResult<index_type, criteria_type>;
		using matelemresult_type_ptr = std::shared_ptr<matelemresult_type>;
		using matelem_type = MatElem<index_type, distance_type, criteria_type>;
		using matord_type = MatOrd<index_type, distance_type, criteria_type>;
		using matordresult_type = std::pair<matelemresult_type_ptr, matelemresult_type_ptr>;
	private:
		std::shared_ptr<matelem_type> m_colmat;
		std::shared_ptr<matelem_type> m_rowmat;
	public:
		MatOrd(cancellableflag_type *pf = nullptr, Backgrounder *pBack = nullptr) :CancellableObject(pf,pBack) {}
		template <typename FLOATTYPE, typename STRINGTYPE>
		MatOrd(MatData<index_type, FLOATTYPE, distance_type, STRINGTYPE> *pMat,
			cancellableflag_type *pf = nullptr, Backgrounder *pBack = nullptr) :
			CancellableObject(pf, pBack) {
			assert(pMat != nullptr);
			assert(pMat->is_valid());
			this->m_colmat = std::make_shared<matelem_type>(pMat, DispositionType::dispCol, pf, pBack);
			assert(this->m_colmat.get() != nullptr);
			this->m_rowmat = std::make_shared<matelem_type>(pMat, DispositionType::dispRow, pf, pBack);
			assert(this->m_rowmat.get() != nullptr);
		}
		MatOrd(const matord_type &other) :CancellableObject(other), m_colmat(other.m_colmat),
			m_rowmat(other.m_rowmat) {}
		matord_type & operator=(const matord_type &other) {
			if (this != other) {
				CancellableObject::operator=(other);
				this->m_colmat = other.m_colmat;
				this->m_rowmat = other.m_rowmat;
			}
			return (*this);
		}
		virtual ~MatOrd() {}
	public:
		bool is_valid(void) const {
			return ((this->m_colmat.get() != nullptr) && this->m_colmat->is_valid() &&
				(this->m_rowmat != nullptr) && this->m_rowmat->is_valid());
		}// is_valid
		matelemresult_type_ptr arrange_cols(void) {
			assert(this->m_colmat.get() != nullptr);
			assert(this->m_colmat->is_valid());
			return (this->m_colmat->arrange());
		}// arrange_cols
		matelemresult_type_ptr arrange_rows(void) {
			assert(this->m_rowmat.get() != nullptr);
			assert(this->m_rowmat->is_valid());
			return (this->m_rowmat->arrange());
		}// arrange_cols
		matordresult_type arrange_all(void) {
			using future_t = std::future<matelemresult_type_ptr>;
			assert(this->is_valid());
			matordresult_type r(std::make_pair(matelemresult_type_ptr(), matelemresult_type_ptr()));
			future_t f = std::async([this]()->matelemresult_type_ptr {
				return (this->arrange_cols());
			});
			r.first = this->arrange_rows();
			r.second = f.get();
			return (r);
		}// arrange_all

		template <typename FLOATTYPE, typename STRINGTYPE>
		static matordresult_type st_arrange_all(MatData<index_type, FLOATTYPE, distance_type, STRINGTYPE> *pMat,
			cancellableflag_type *pf = nullptr, Backgrounder *pBack = nullptr) {
			matord_type oMat(pMat, pf, pBack);
			return oMat.arrange_all();
		}// st_arrange_all
		template <typename FLOATTYPE, typename STRINGTYPE>
		static  matordresult_type st_arrange_all_hierar(MatData<index_type, FLOATTYPE, distance_type, STRINGTYPE> *pMat,
			cancellableflag_type *pf = nullptr, Backgrounder *pBack = nullptr) {
			using future_t = std::future<matelemresult_type_ptr>;
			matordresult_type oRet = std::make_pair(matelemresult_type_ptr(), matelemresult_type_ptr());
			future_t f = std::async([pMat, pf, pBack, &oRet]()->matelemresult_type_ptr {
				MatTree<index_type, FLOATTYPE, distance_type, STRINGTYPE> oTree(pMat, DispositionType::dispCol, pf);
				if (!oTree.aggreg()) {
					return (matelemresult_type_ptr());
				}
				sizets_vector oinds;
				oTree.get_ids(oinds);
				const size_t n = oinds.size();
				for (size_t i = 0; i < n; ++i) {
					index_type ix = oinds[i];
					assert(ix > 0);
					oinds[i] = (index_type)(ix - 1);
				}// i
				matelem_type oMat(pMat, DispositionType::dispCol, pf, pBack);
				criteria_type c = oMat.criteria();
				oMat.set(c, oinds);
				return (oMat.arrange());
			});
			{
				MatTree<index_type, FLOATTYPE, distance_type, STRINGTYPE> oTree(pMat, DispositionType::dispRow, pf);
				if (!oTree.aggreg()) {
					return (oRet);
				}
				sizets_vector oinds;
				oTree.get_ids(oinds);
				const size_t n = oinds.size();
				for (size_t i = 0; i < n; ++i) {
					index_type ix = oinds[i];
					assert(ix > 0);
					oinds[i] = (index_type)(ix - 1);
				}// i
				matelem_type oMat(pMat, DispositionType::dispRow, pf, pBack);
				criteria_type c = oMat.criteria();
				oMat.set(c, oinds);
				oRet.first = oMat.arrange();
			}// rows
			oRet.second = f.get();
			return (oRet);
		}// create
	}; // class MatOrd
	   //////////////////////////////////
}// namespace info
 ///////////////////////////////
#endif // __MATELEM_H__
