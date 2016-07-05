#pragma once
#ifndef __DISTANCE_MAP_H__
#define __DISTANCE_MAP_H__
/////////////////////////////////
#include "gendefs.h"
///////////////////////////////////////
namespace info {
	////////////////////////////////////
	template<typename INDEXTYPE, typename DISTANCETYPE>
	class DistanceMap {
	public:
		using index_type = INDEXTYPE;
		using distance_type = DISTANCETYPE;
		//
		using first_map_type = std::map<index_type, distance_type>;
		using map_type = std::map<index_type, first_map_type>;
		using set_type = std::set<index_type>;
		using distancemap_type = DistanceMap<index_type, distance_type>;
	private:
		map_type m_map;
		set_type m_set;
	public:
		DistanceMap() {}
		DistanceMap(const distancemap_type &other) :
			m_map(other.m_map), m_set(other.m_set) {
		}
		distancemap_type & operator=(const distancemap_type &other) {
			if (this != &other) {
				this->m_map = other.m_map;
				this->m_set = other.m_set;
			}
			return (*this);
		}
		virtual ~DistanceMap() {
		}
	public:
		bool is_valid(void) const {
			return ((this->m_map.size() > 0) && (this->m_set.size() > 0));
		}// i_valid
		template <typename XU, typename CRITERIATYPE>
		void compute_criteria(const std::vector<XU> &oInds, CRITERIATYPE &res) const {
			double ss = 0;
			const size_t n = oInds.size();
			for (size_t i = 1; i < n; ++i) {
				index_type ii = (XU)oInds[i - 1];
				index_type jj = (XU)oInds[i];
				distance_type val = 0;
				if (this->get(ii, jj, val)) {
					ss += (double)val;
				}
			}// i
			res = (CRITERIATYPE)ss;
		}// compute_criteria
		template <typename XU>
		void get_indexes(std::vector<XU> &v) const {
			v.clear();
			const set_type &oSet = this->m_set;
			for (auto it = oSet.begin(); it != oSet.end(); ++it) {
				XU val = (XU)(*it);
				v.push_back(val);
			} // it
		} // get_indexes
		void clear(void) {
			this->m_map.clear();
			this->m_set.clear();
		} // clear
		template <typename XU>
		bool has_entry(const XU ii, const XU jj) const {
			const set_type &oSet = this->m_set;
			index_type i = (XU)ii;
			index_type j = (XU)jj;
			return ((oSet.find(i) != oSet.end()) && (oSet.find(j) != oSet.end()));
		} // has_entry
		template <typename XU, typename XF>
		bool get(const XU ii, const XU jj, XF &val) const {
			index_type i = (index_type)ii, j = (index_type)jj;
			if (i > j) {
				index_type t = i;
				i = j;
				j = t;
			}
			const map_type &oMap = this->m_map;
			auto it = oMap.find(i);
			if (it == oMap.end()) {
				return (false);
			}
			const first_map_type & m = (*it).second;
			auto jt = m.find(j);
			if (jt != m.end()) {
				val = (XF) ((*jt).second);
				return (true);
			}
			return (false);
		} // get
		template <typename XU, typename XF>
		bool add(const XU ii, const XU jj, const XF val) {
			if (ii == jj) {
				return (false);
			}
			index_type i = (index_type)ii, j = (index_type)jj;
			if (i > j) {
				index_type t = i;
				i = j;
				j = t;
			}
			distance_type zval = (distance_type)val;
			map_type &oMap = this->m_map;
			set_type &oSet = this->m_set;
			auto it = oMap.find(i);
			if (it == oMap.end()) {
				first_map_type m;
				m[j] = zval;
				oMap[i] = m;
				oSet.insert(i);
				oSet.insert(j);
			}
			else {
				first_map_type &m = (*it).second;
				m[j] = zval;
				oSet.insert(j);
			}
			return (true);
		} // add
		template<typename XU, typename X>
		bool recode(DistanceMap<XU, X> &oDest, X dMax = 1000, X dMin = 0) const {
			if (dMax <= dMin) {
				return (false);
			}
			distance_type vMin = 0, vMax = 0;
			bool bFirst = true;
			const map_type &oMap = this->m_map;
			for (auto it = oMap.begin(); it != oMap.end(); ++it) {
				const first_map_type &m = (*it).second;
				for (auto jt = m.begin(); jt != m.end(); ++jt) {
					distance_type v = (*jt).second;
					if (bFirst) {
						vMin = v;
						vMax = v;
						bFirst = false;
					}
					else if (v < vMin) {
						vMin = v;
					}
					else if (v > vMax) {
						vMax = v;
					}
				} // jt
			} // it
			if (vMax <= vMin) {
				return (false);
			}
			double delta = (double)(dMax - dMin) / (double)(vMax - vMin);
			oDest.clear();
			for (auto it = oMap.begin(); it != oMap.end(); ++it) {
				index_type ii = (*it).first;
				const first_map_type &m = (*it).second;
				for (auto jt = m.begin(); jt != m.end(); ++jt) {
					index_type jj = (*jt).first;
					distance_type v = (*jt).second;
					double x = (((double)v - vMin) * delta) + (double)dMin;
					X val = (X)x;
					oDest.add(ii, jj, val);
				} // jt
			} // it
			return (true);
		} // recode
		template<typename XU, typename X>
		std::future<bool> recodeAsync(DistanceMap<XU, X> &oDest, X dMax = 1000, X dMin = 0) const {
			return std::async([this, dMax, dMin, &oDest]()->bool {
				return (this->recode(oDest, dMax, dMin));
			});
		}// recodeAsync
	};
	// class DistanceMap<T>
	/////////////////////////////////////
}// namespace info
///////////////////////////////////
#endif // __DISTANCE_MAP_H__

