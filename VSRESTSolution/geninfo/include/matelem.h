#pragma once
#ifndef __MATELEM_H__
#define __MATELEM_H__
/////////////////////////////
#include "matdata.h"
//////////////////////
#include <atomic>
#include <queue>
////////////////////////
namespace info {
	////////////////////////////////
	template <typename DISTANCETYPE, typename CRITERIATYPE>
	class MatElem {
	public:
		using criteria_type = CRITERIATYPE;
		using sizets_vector = std::vector<size_t>;
		using DistanceMapType = DistanceMap<size_t, DISTANCETYPE>;
		using pair_type = std::pair<size_t, size_t>;
		using pairs_type_vector = std::queue<pair_type>;
		using index_ptr_type = std::shared_ptr<sizets_vector>;
		using result_type = std::pair<criteria_type, index_ptr_type>;
		using task_type = pplx::task<result_type>;
		using MatElemType = MatElem<DISTANCETYPE, CRITERIATYPE>;
	public:
		MatElem(const DistanceMapType *pMap) :m_pdist(pMap), m_cancel(false), m_crit(0) {
			assert(this->m_pdist != nullptr);
			pMap->get_indexes(this->m_indexes);
			pMap->compute_criteria(this->m_indexes, this->m_crit);
		}// MatElem
		virtual ~MatElem() {}
	public:
		result_type get_result(void) {
			index_ptr_type oPtr = std::make_shared<sizets_vector>(this->m_indexes);
			return std::make_pair(this->m_crit, oPtr);
		}// get_result
		void cancel(void) {
			this->m_cancel.store(true);
		}
		result_type arrange(void) {
			bool done = false;
			while (!done) {
				if (this->m_cancel.load()) {
					done = true;
					break;
				}
				pairs_type_vector v;
				criteria_type bestCrit = this->m_crit;
				if (!this->find_best_permutations(v, bestCrit)) {
					done = true;
					break;
				}
				if (v.empty()) {
					done = true;
					break;
				}
				pair_type p = v.front();
				v.pop();
				size_t iFirst = p.first;
				size_t jFirst = p.second;
				std::vector<task_type> mtasks;
				while (!v.empty()) {
					pair_type px = v.front();
					size_t i = px.first;
					size_t j = px.second;
					v.pop();
					task_type ts([this, i, j, bestCrit]()->result_type {
						MatElemType oMat(*this);
						size_t t = oMat.m_indexes[i];
						oMat.m_indexes[i] = oMat.m_indexes[j];
						oMat.m_indexes[j] = t;
						oMat.m_crit = bestCrit;
						return (oMat.arrange());
					});
					mtasks.push_back(ts);
				}// not empty
				sizets_vector &xind = this->m_indexes;
				size_t t = xind[iFirst];
				xind[iFirst] = xind[jFirst];
				xind[jFirst] = t;
				this->m_crit = bestCrit;
				for (auto it = mtasks.begin(); it != mtasks.end(); ++it) {
					task_type ts = *it;
					result_type pp = ts.get();
					const criteria_type crit = pp.first;
					if (crit < this->m_crit) {
						index_ptr_type oind = pp.second;
						const sizets_vector *pv = oind.get();
						if (pv != nullptr) {
							this->m_crit = crit;
							this->m_indexes = *pv;
						}// pv
					}
				}// it
			}// not done
			return (this->get_result());
		}// arrange
		task_type arrangeAsync(void) {
			return pplx::task<result_type>([this]()->result_type {
				return (this->arrange());
			});
		}// arrangeAsync
	protected:
		MatElem(const MatElemType &other) :m_pdist(other.m_pdist), m_cancel(false), m_crit(other.m_crit),
			m_indexes(other.m_indexes) {}
		MatElemType & operator=(const MatElemType &other) {
			if (this != &other) {
				this->m_pdist = other.m_pdist;
				this->m_crit = other.m_crit;
				this->m_indexes = other.m_indexes;
			}
			return (*this);
		}
	private:
		bool find_best_permutations(pairs_type_vector &v, criteria_type &bestCrit) const {
			const sizets_vector &inds = this->m_indexes;
			const size_t n = inds.size();
			const DistanceMapType *pMap = this->m_pdist;
			for (size_t i = 0; i < n; ++i) {
				if (this->m_cancel.load()) {
					return (false);
				}
				for (size_t j = 0; j < i; ++j) {
					if (this->m_cancel.load()) {
						return (false);
					}
					sizets_vector xind(inds);
					size_t t = xind[i];
					xind[i] = xind[j];
					xind[j] = t;
					criteria_type c = 0;
					pMap->compute_criteria(xind, c);
					if (c < bestCrit) {
						if (!v.empty()) {
							while (!v.empty()) {
								v.pop();
							}
						}
						v.push(std::make_pair(j, i));
						bestCrit = c;
					}
					else if ((c == bestCrit) && (!v.empty())) {
						v.push(std::make_pair(j, i));
					}
				}// j
			}// i
			if (this->m_cancel.load()) {
				return (false);
			}
			return (!v.empty());
		}//find_best_permutations
	private:
		const DistanceMapType *m_pdist;
		std::atomic<bool> m_cancel;
		criteria_type m_crit;
		sizets_vector m_indexes;
	}; // class MatElem
	//////////////////////////////////
}// namespace info
/////////////////////////////
#endif // __MATELEM_H__
