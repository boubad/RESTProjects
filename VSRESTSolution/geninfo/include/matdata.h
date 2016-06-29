/*
 * matdata.h
 *
 *  Created on: 27 juin 2016
 *      Author: boubad
 */

#pragma once
#ifndef __MATDATA_H__
#define __MATDATA_H__
 ///////////////////////////////////////
#include <cpprest/details/basic_types.h>
#include <pplx/pplxtasks.h>
///////////////////////////////////////
#include <cassert>
#include <memory>
#include <vector>
#include <map>
#include <set>
//////////////////////////////////////////////////
namespace info {
	////////////////////////////////////////
	template<typename U, typename T>
	class DistanceMap {
	public:
		using ints_vector = std::vector<U>;
		using first_map_type = std::map<U, T>;
		using map_type = std::map<U, first_map_type>;
		using set_type = std::set<U>;
		using DistanceMapType = DistanceMap<U, T>;
	private:
		map_type m_map;
		set_type m_set;
	public:
		DistanceMap() {
		}
		DistanceMap(const DistanceMapType &other) :
			m_map(other.m_map), m_set(other.m_set) {
		}
		DistanceMapType & operator=(const DistanceMapType &other) {
			if (this != &other) {
				this->m_map = other.m_map;
				this->m_set = other.m_set;
			}
			return (*this);
		}
		virtual ~DistanceMap() {
		}
	public:
		template <typename X, typename F>
		void compute_criteria(const std::vector<X> &oInds, F &res) const {
			double ss = 0;
			size_t nc = 0;
			const size_t n = oInds.size();
			for (size_t i = 1; i < n; ++i) {
				U ii = (U)oInds[i - 1];
				U jj = (U)oInds[i];
				T val;
				if (this->get(ii, jj, val)) {
					ss += (double)val;
					++nc;
				}
			}// i
			if (nc > 1) {
				ss /= nc;
			}
			res = (F)ss;
		}// compute_criteria
		template <typename X>
		void get_indexes(std::vector<X> &v) const {
			v.clear();
			const set_type &oSet = this->m_set;
			for (auto it = oSet.begin(); it != oSet.end(); ++it) {
				X val = (X)(*it);
				v.push_back(val);
			} // it
		} // get_indexes
		void clear(void) {
			this->m_map.clear();
			this->m_set.clear();
		} // clear
		template <typename X>
		bool has_entry(const X ii, const X jj) const {
			const set_type &oSet = this->m_set;
			U i = (U)ii;
			U j = (U)jj;
			return ((oSet.find(i) != oSet.end()) && (oSet.find(j) != oSet.end()));
		} // has_entry
		template <typename X>
		bool get(const X ii, const X jj, T &val) const {
			U i = (U)ii, j = (U)jj;
			if (i > j) {
				U t = i;
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
				val = (*jt).second;
				return (true);
			}
			return (false);
		} // get
		template <typename X>
		bool add(const X ii, const X jj, T val) {
			if (ii == jj) {
				return (false);
			}
			U i = (U)ii, j = (U)jj;
			if (i > j) {
				U t = i;
				i = j;
				j = t;
			}
			map_type &oMap = this->m_map;
			set_type &oSet = this->m_set;
			auto it = oMap.find(i);
			if (it == oMap.end()) {
				first_map_type m;
				m[j] = val;
				oMap[i] = m;
				oSet.insert(i);
				oSet.insert(j);
			}
			else {
				first_map_type &m = (*it).second;
				m[j] = val;
				oSet.insert(j);
			}
			return (true);
		} // add
		template<typename XU, typename X>
		bool recode(DistanceMap<XU, X> &oDest, X dMax = 1000, X dMin = 0) const {
			if (dMax <= dMin) {
				return (false);
			}
			T vMin = 0, vMax = 0;
			bool bFirst = true;
			const map_type &oMap = this->m_map;
			for (auto it = oMap.begin(); it != oMap.end(); ++it) {
				const first_map_type &m = (*it).second;
				for (auto jt = m.begin(); jt != m.end(); ++jt) {
					T v = (*jt).second;
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
				U ii = (*it).first;
				const first_map_type &m = (*it).second;
				for (auto jt = m.begin(); jt != m.end(); ++jt) {
					U jj = (*jt).first;
					T v = (*jt).second;
					double x = (((double)v - vMin) * delta) + (double)dMin;
					X val = (X)x;
					oDest.add(ii, jj, val);
				} // jt
			} // it
			return (true);
		} // recode
		template<typename XU, typename X>
		pplx::task<bool> recodeAsync(DistanceMap<XU, X> &oDest, X dMax = 1000, X dMin = 0) const {
			return pplx::task<bool>([this, dMax, dMin, &oDest]()->bool {
				return (this->recode(oDest, dMax, dMin));
			});
		}// recodeAsync
	};
	// class DistanceMap<T>
	/////////////////////////////////////
	template <typename T>
	class MatArgs {
	public:
		using data_type = T;
		using string_type = utility::string_t;
		using strings_vector = std::vector<string_type>;
		using data_vector = std::vector<data_type>;
	private:
		size_t m_nrows;
		size_t m_ncols;
		data_vector m_data;
		strings_vector m_rownames;
		strings_vector m_colnames;
	protected:
		bool set_names(const strings_vector *pRowNames = nullptr,
			const strings_vector *pColNames = nullptr) {
			const size_t nRows = this->m_nrows;
			const size_t nCols = this->m_ncols;
			if ((nRows < 1) || (nCols < 1)) {
				return (false);
			}
			strings_vector rownames(nRows), colnames(nCols);
			if ((pRowNames != nullptr) && (pRowNames->size() >= nRows)) {
				const strings_vector &v = *pRowNames;
				for (size_t i = 0; i < nRows; ++i) {
					string_type s = v[i];
					rownames[i] = s;
				}	// i
			}
			else {
				for (size_t i = 0; i < nRows; ++i) {
					utility::stringstream_t os;
					os << U("ind") << (i + 1);
					rownames[i] = os.str();
				}	// i
			}
			if ((pColNames != nullptr) && (pColNames->size() >= nCols)) {
				const strings_vector &v = *pColNames;
				for (size_t i = 0; i < nCols; ++i) {
					string_type s = v[i];
					colnames[i] = s;
				}	// i
			}
			else {
				for (size_t i = 0; i < nCols; ++i) {
					utility::stringstream_t os;
					os << U("var") << (i + 1);
					colnames[i] = os.str();
				}	// i
			}
			this->m_rownames = rownames;
			this->m_colnames = colnames;
			return (true);
		}	// setnames
	public:
		template <typename X>
		MatArgs(size_t nRows, size_t nCols, const X &xdata, const strings_vector *pRowNames = nullptr,
			const strings_vector *pColNames = nullptr) :m_nrows(nRows), m_ncols(nCols) {
			const size_t nn = (size_t)(nCols * nRows);
			this->m_data.resize(nn);
			for (size_t i = 0; i < nn; ++i) {
				(this->m_data)[i] = (T)xdata[i];
			}// i
			this->set_names(pRowNames, pColNames);
		}// MatArgs
		MatArgs(const MatArgs<T> &other) :m_nrows(other.m_nrows), m_ncols(other.m_ncols),
			m_data(other.m_data), m_rownames(other.m_rownames), m_colnames(other.m_colnames) {}
		MatArgs<T> & operator=(const MatArgs<T> &other) {
			if (this != &other) {
				this->m_nrows = other.m_nrows;
				this->m_ncols = other.m_ncols;
				this->m_data = other.m_data;
				this->m_rownames = other.m_rownames;
				this->m_colnames = other.m_colnames;
			}
			return (*this);
		}
		virtual ~MatArgs() {}
	public:
		bool is_valid(void) const {
			return ((this->m_nrows > 0) && (this->m_ncols > 0));
		}	// is_valid
		size_t rows(void) const {
			return (this->m_nrows);
		}	// rows
		size_t cols(void) const {
			return (this->m_ncols);
		}	// cols
		const strings_vector &row_names(void) const {
			return (this->m_rownames);
		}
		const strings_vector &col_names(void) const {
			return (this->m_colnames);
		}
		data_type get_value(const size_t irow, const size_t icol) const {
			assert(irow < this->m_nrows);
			const size_t nCols = this->m_ncols;
			assert(icol < nCols);
			return ((this->m_data)[irow * nCols + icol]);
		}	// get_value
		data_type operator()(const size_t irow, const size_t icol) const {
			return (this->get_value(irow, icol));
		}	// operator()
	}; // class MatArgs<T>
	///////////////////////////////////
	template <typename INDEXTYPE, typename FLOATTYPE, typename DISTANCETYPE>
	class MatData {
	public:
		using index_type = INDEXTYPE;
		using data_type = FLOATTYPE;
		using distance_type = DISTANCETYPE;
		using string_type = utility::string_t;
		//
		using strings_vector = std::vector<string_type>;
		using doubles_vector = std::vector<data_type>;
		using DistanceMapType = DistanceMap<index_type, distance_type>;
		using PDistanceMapType = DistanceMapType *;
		using MatDataType = MatData<index_type, data_type, distance_type>;
		using MatDataPtr = std::shared_ptr<MatDataType>;
	private:
		size_t m_nrows;
		size_t m_ncols;
		std::shared_ptr<data_type> m_data;
		strings_vector m_rownames;
		strings_vector m_colnames;
		std::shared_ptr<DistanceMapType> m_rowdist;
		std::shared_ptr<DistanceMapType> m_coldist;
	public:
		MatData() :
			m_nrows(0), m_ncols(0) {
		}
		MatData(const MatDataType &other) :
			m_nrows(other.m_nrows), m_ncols(other.m_ncols), m_data(
				other.m_data), m_rownames(other.m_rownames), m_colnames(
					other.m_colnames), m_rowdist(other.m_rowdist), m_coldist(
						other.m_coldist) {
		}
		MatDataType & operator=(const MatDataType &other) {
			if (this != &other) {
				this->m_nrows = other.m_nrows;
				this->m_ncols = other.m_ncols;
				this->m_data = other.m_data;
				this->m_rownames = other.m_rownames;
				this->m_colnames = other.m_colnames;
				this->m_rowdist = other.m_rowdist;
				this->m_coldist = other.m_coldist;
			}
			return (*this);
		}
		virtual ~MatData() {
		}
	public:
		bool is_valid(void) const {
			return ((this->m_nrows > 0) && (this->m_ncols > 0));
		}	// is_valid
		size_t rows(void) const {
			return (this->m_nrows);
		}	// rows
		size_t cols(void) const {
			return (this->m_ncols);
		}	// cols
		data_type get_value(const size_t irow, const size_t icol) const {
			assert(irow < this->m_nrows);
			const size_t nCols = this->m_ncols;
			assert(icol < nCols);
			const data_type *p = this->m_data.get();
			assert(p != nullptr);
			return (p[irow * nCols + icol]);
		}	// get_value
		data_type operator()(const size_t irow, const size_t icol) const {
			return (this->get_value(irow, icol));
		}	// operator()
		string_type row_name(const size_t irow) const {
			const strings_vector &v = this->m_rownames;
			assert((irow < this->m_nrows) && (irow < v.size()));
			return (v[irow]);
		}	// row_name
		string_type col_name(const size_t icol) const {
			const strings_vector &v = this->m_colnames;
			assert((icol < this->m_ncols) && (icol < v.size()));
			return (v[icol]);
		}	// row_name
		bool set_names(const strings_vector *pRowNames = nullptr,
			const strings_vector *pColNames = nullptr) {
			const size_t nRows = this->m_nrows;
			const size_t nCols = this->m_ncols;
			if ((nRows < 1) || (nCols < 1)) {
				return (false);
			}
			strings_vector rownames(nRows), colnames(nCols);
			if ((pRowNames != nullptr) && (pRowNames->size() >= nRows)) {
				const strings_vector &v = *pRowNames;
				for (size_t i = 0; i < nRows; ++i) {
					string_type s = v[i];
					rownames[i] = s;
				}	// i
			}
			else {
				for (size_t i = 0; i < nRows; ++i) {
					utility::stringstream_t os;
					os << U("ind") << (i + 1);
					rownames[i] = os.str();
				}	// i
			}
			if ((pColNames != nullptr) && (pColNames->size() >= nCols)) {
				const strings_vector &v = *pColNames;
				for (size_t i = 0; i < nCols; ++i) {
					string_type s = v[i];
					colnames[i] = s;
				}	// i
			}
			else {
				for (size_t i = 0; i < nCols; ++i) {
					utility::stringstream_t os;
					os << U("var") << (i + 1);
					colnames[i] = os.str();
				}	// i
			}
			this->m_rownames = rownames;
			this->m_colnames = colnames;
			return (true);
		}	// setnames
	public:
		template<typename X>
		bool init(std::shared_ptr<MatArgs<X> > oArgs) {
			const MatArgs<X> *pArgs = oArgs.get();
			if (pArgs == nullptr) {
				return (false);
			}
			if (!pArgs->is_valid()) {
				return (false);
			}
			const size_t nRows = pArgs->rows();
			const size_t nCols = pArgs->cols();
			bool bRet = false;
			try {
				const size_t nn = (size_t)(nRows * nCols);
				std::shared_ptr<data_type> ff(new data_type[nn]);
				data_type *pf = ff.get();
				if (pf == nullptr) {
					return (bRet);
				}
				for (size_t i = 0; i < nCols; ++i) {
					double vmin = 0, vmax = 0;
					for (size_t j = 0; j < nRows; ++j) {
						double x = (double)pArgs->get_value(j, i);
						if (j == 0) {
							vmin = x;
							vmax = x;
						}
						else if (x < vmin) {
							vmin = x;
						}
						else if (x > vmax) {
							vmax = x;
						}
					}	// j
					if (vmin >= vmax) {
						return (bRet);
					}
					const double delta = vmax - vmin;
					for (size_t j = 0; j < nRows; ++j) {
						pf[j * nCols + i] = (data_type)((pArgs->get_value(j, i) - vmin) / delta);
					}	// j
				}	// i
				this->m_nrows = nRows;
				this->m_ncols = nCols;
				this->m_data = ff;
				this->m_rownames = pArgs->row_names();
				this->m_colnames = pArgs->col_names();
				this->m_rowdist.reset();
				this->m_coldist.reset();
				bRet = true;
			}
			catch (...) {}
			return (bRet);
		}	// init
		PDistanceMapType get_rows_distances_map(void) {
			PDistanceMapType pRet = this->m_rowdist.get();
			if (pRet == nullptr) {
				this->getRowDist();
				pRet = this->m_rowdist.get();
			}
			return (pRet);
		}							// get_rows_distances_map
		PDistanceMapType get_cols_distances_map(void) {
			PDistanceMapType pRet = this->m_coldist.get();
			if (pRet == nullptr) {
				this->getColDist();
				pRet = this->m_coldist.get();
			}
			return (pRet);
		}							 // get_cols_distances_map
		template<typename X>
		static MatDataPtr create(std::shared_ptr<MatArgs<X> > oArgs) {
			MatDataPtr oRet = std::make_shared<MatDataType>();
			MatDataType *pRet = oRet.get();
			if (pRet != nullptr) {
				pRet->init(oArgs);
				pRet->getRowDist();
				pRet->getColDist();
			}// pRet
			return (oRet);
		}		// create
				/////////////////////////////////////
	protected:
		void getRowDist(void) {
			const size_t nRows = this->m_nrows;
			const size_t nCols = this->m_ncols;
			const FLOATTYPE *pf = this->m_data.get();
			if ((nRows > 0) && (nCols > 0) && (pf != nullptr)) {
				DistanceMap<size_t, double> oDist;
				for (size_t i = 0; i < nRows; ++i) {
					const FLOATTYPE *p1 = pf + i * nCols;
					for (size_t j = 0; j < i; ++j) {
						const FLOATTYPE *p2 = pf + (j * nCols);
						double s = 0;
						for (size_t k = 0; k < nCols; ++k) {
							double x = p1[k] - p2[k];
							s += x * x;
						}	// k
						double val = s / nCols;
						oDist.add(j, i, val);
					}	// j
				}	// i
					//
				this->m_rowdist = std::make_shared<DistanceMapType>();
				DistanceMapType *pRet = this->m_rowdist.get();
				if (pRet != nullptr) {
					oDist.recode(*pRet);
				}
			}
		}	// getRowDistAsync
		void getColDist(void) {
			const size_t nRows = this->m_nrows;
			const size_t nCols = this->m_ncols;
			const FLOATTYPE *pf = this->m_data.get();
			if ((nRows > 0) && (nCols > 0) && (pf != nullptr)) {
				DistanceMap<size_t, double> oDist;
				for (size_t i = 0; i < nCols; ++i) {
					for (size_t j = 0; j < i; ++j) {
						double s = 0;
						for (size_t k = 0; k < nRows; ++k) {
							double x = pf[k * nCols + i] - pf[k * nCols + j];
							s += x * x;
						}							// k
						double val = s / nRows;
						oDist.add(j, i, val);
					}							// j
				}					//i							// i
				this->m_coldist = std::make_shared<DistanceMapType>();
				DistanceMapType *pRet = this->m_coldist.get();
				if (pRet != nullptr) {
					oDist.recode(*pRet);
				}
			}							 // ok
		}							 // getColistAsync
	};
	// class MatData
	///////////////////////////////////////
}// namespace info
////////////////////////////////////////////
#endif // !__MATDATA_H__
