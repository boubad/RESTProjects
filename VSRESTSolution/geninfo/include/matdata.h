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
	template <typename T>
	class DistanceMap {
	public:
		using sizets_vector = std::vector<size_t>;
		using first_map_type = std::map<size_t, T>;
		using map_type = std::map<size_t, first_map_type>;
		using set_type = std::set<T>;
	private:
		map_type m_map;
		set_type m_set;
	public:
		DistanceMap() {}
		DistanceMap(const DistanceMap<T> &other) :m_map(other.m_map), m_set(other.m_set) {}
		DistanceMap<T> & operator=(const DistanceMap<T> &other) {
			if (this != &other) {
				this->m_map = other.m_map;
				this->m_set = other.m_set;
			}
			return (*this);
		}
		virtual ~DistanceMap() {}
	public:
		void get_indexes(sizets_vector &v) const {
			v.clear();
			const set_type &oSet = this->m_set;
			for (auto it = oSet.begin(); it != oSet.end(); ++it) {
				v.push_back(*it);
			}// it
		}// get_indexes
		void clear(void) {
			this->m_map.clear();
			this->m_set.clear();
		}// clear
		bool has_entry(const size_t ii, const size_t jj) const {
			const set_type &oSet = this->m_set;
			return ((oSet.find(ii) != oSet.end()) && (oSet.find(jj) != oSet.end()));
		}// has_entry
		bool get(const size_t ii, const size_t jj, T &val) const {
			size_t i = ii, j = jj;
			if (i > j) {
				size_t t = i;
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
		}// get
		bool add(const size_t ii, const size_t jj, T val) {
			if (ii == jj) {
				return (false);
			}
			size_t i = ii, j = jj;
			if (i > j) {
				size_t t = i;
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
		}// add
		template <typename X>
		bool recode(DistanceMap<X> &oDest, X dMax = 1000, X dMin = 0) const {
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
				}// jt
			}// it
			if (vMax <= vMin) {
				return (false);
			}
			double delts = (double)(dMax - dMin) / (double)(vMax - vMin);
			oDest.clear();
			for (auto it = oMap.begin(); it != oMap.end(); ++it) {
				size_t ii = (*it).first;
				const first_map_type &m = (*it).second;
				for (auto jt = m.begin(); jt != m.end(); ++jt) {
					size_t jj = (*jt).first;
					T v = (*jt).second;
					double x = (((double)v - vMin) * delta) + (double)dMin;
					X val = (X)x;
					oDest.add(ii, jj, val);
				}// jt
			}// it
			return (true);
		}// recode
	};// class DistanceMap<T>
	///////////////////////////////////
	class MatData {
	public:
		using string_type = utility::string_t;
		using strings_vector = std::vector<string_type>;
		using doubles_vector = std::vector<double>;
		using DistanceMapType = DistanceMap<double>;
		using PDistanceMapType = DistanceMapType *;
		using MatDataPtr = std::shared_ptr<MatData>;
	private:
		size_t m_nrows;
		size_t m_ncols;
		std::shared_ptr<double> m_data;
		std::shared_ptr<strings_vector>  m_rownames;
		std::shared_ptr<strings_vector> m_colnames;
		std::shared_ptr<DistanceMapType> m_rowdist;
		std::shared_ptr<DistanceMapType> m_coldist;
	public:
		MatData() : m_nrows(0), m_ncols(0) {}
		MatData(const MatData &other) :m_nrows(other.m_nrows), m_ncols(other.m_ncols),
			m_data(other.m_data), m_rownames(other.m_rownames), m_colnames(other.m_colnames),
			m_rowdist(other.m_rowdist), m_coldist(other.m_coldist) {}
		MatData & operator=(const MatData &other) {
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
		virtual ~MatData() {}
	public:
		bool is_valid(void) const {
			return ((this->m_nrows > 0) && (this->m_ncols > 0));
		}// is_valid
		size_t rows(void) const {
			return (this->m_nrows);
		}// rows
		size_t cols(void) const {
			return (this->m_ncols);
		}// cols
		double get_value(const size_t irow, const size_t icol) const {
			assert(irow < this->m_nrows);
			const size_t nCols = this->m_ncols;
			assert(icol < nCols);
			const double *p = this->m_data.get();
			assert(p != nullptr);
			return (p[irow * nCols + icol]);
		}// get_value
		double operator()(const size_t irow, const size_t icol) const {
			return (this->get_value(irow, icol));
		}// operator()
		string_type row_name(const size_t irow) const {
			const strings_vector *pv = this->m_rownames.get();
			assert(pv != nullptr);
			const strings_vector &v = *pv;
			assert((irow < this->m_nrows) && (irow < v.size()));
			return (v[irow]);
		}// row_name
		string_type col_name(const size_t icol) const {
			const strings_vector *pv = this->m_colnames.get();
			assert(pv != nullptr);
			const strings_vector &v = *pv;
			assert((icol < this->m_ncols) && (icol < v.size()));
			return (v[icol]);
		}// row_name
		bool set_names(const strings_vector *pRowNames = nullptr, const strings_vector *pColNames = nullptr) {
			const size_t nRows = this->m_nrows;
			const size_t nCols = this->m_ncols;
			if ((nRows < 1) || (nCols < 1)) {
				return (false);
			}
			std::shared_ptr<strings_vector> rownames = std::make_shared<strings_vector>();
			std::shared_ptr<strings_vector> colnames = std::make_shared<strings_vector>();
			strings_vector *pr = rownames.get();
			strings_vector *pc = colnames.get();
			if ((pr == nullptr) || (pc == nullptr)) {
				return (false);
			}
			strings_vector &nr = *pr;
			nr.resize(nRows);
			strings_vector &nc = *pc;
			nc.resize(nCols);
			if ((pRowNames != nullptr) && (pRowNames->size() >= nRows)) {
				const strings_vector &v = *pRowNames;
				for (size_t i = 0; i < nRows; ++i) {
					string_type s = v[i];
					nr[i] = s;
				}// i
			}
			else {
				for (size_t i = 0; i < nRows; ++i) {
					utility::stringstream_t os;
					os << U("ind") << (i + 1);
					nr[i] = os.str();
				}// i
			}
			if ((pColNames != nullptr) && (pColNames->size() >= nCols)) {
				const strings_vector &v = *pColNames;
				for (size_t i = 0; i < nCols; ++i) {
					string_type s = (*pColNames)[i];
					nc[i] = s;
				}// i
			}
			else {
				for (size_t i = 0; i < nCols; ++i) {
					utility::stringstream_t os;
					os << U("var") << (i + 1);
					nc[i] = os.str();
				}// i
			}
			this->m_rownames = rownames;
			this->m_colnames = colnames;
			return (true);
		}// setnames
	public:
		template <typename X>
		pplx::task<bool> initAsync(const size_t nRows, const size_t nCols, const X &pxdata) {
			return pplx::task<bool>([this, nRows, nCols, pxdata]()->bool {
				bool bRet = false;
				if ((nRows < 1) || (nCols < 1)) {
					return (bRet);
				}
				try {
					const size_t nn = (size_t)(nRows * nCols);
					std::shared_ptr<double> ff(new double[nn]);
					double *pf = ff.get();
					if (pf == nullptr) {
						return (bRet);
					}
					for (size_t i = 0; i < nCols; ++i) {
						double vmin = 0, vmax = 0;
						for (size_t j = 0; j < nRows; ++j) {
							double x = (double)pxdata[j * nCols + i];
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
						}// j
						if (vmin >= vmax) {
							return (bRet);
						}
						const double delta = vmax - vmin;
						for (size_t j = 0; j < nRows; ++j) {
							const size_t pos = j * nCols + i;
							pf[pos] = (pxdata[pos] - vmin) / delta;
						}// j
					}// i
					this->m_nrows = nRows;
					this->m_ncols = nCols;
					this->m_data = ff;
					this->m_rowdist.reset();
					this->m_coldist.reset();
					bRet = true;
				}
				catch (...) {}
				return (bRet);
			});
		}// iniAsync
		template <typename X>
		bool init(const size_t nRows, const size_t nCols, const X &pxdata) {
			pplx::task<bool> tsk = this->initAsync(nRows, nCols, pxdata);
			bool bRet = tsk.get();
			return (bRet);
		}// init
		pplx::task<PDistanceMapType> getRowDistAsync(void) {
			return pplx::task<PDistanceMapType>([this]()->PDistanceMapType {
				PDistanceMapType pRet = this->m_rowdist.get();
				if (pRet == nullptr) {
					const size_t nRows = this->m_nrows;
					const size_t nCols = this->m_ncols;
					const double *pf = this->m_data.get();
					if ((nRows > 0) && (nCols > 0) && (pf != nullptr)) {
						std::shared_ptr<DistanceMapType> oDist(new DistanceMapType());
						PDistanceMapType pDest = oDist.get();
						if (pDest != nullptr) {
							//
							for (size_t i = 0; i < nRows; ++i) {
								const double *p1 = pf + i * nCols;
								for (size_t j = 0; j < i; ++j) {
									const double *p2 = pf + (j * nCols);
									double s = 0;
									for (size_t k = 0; k < nCols; ++k) {
										double x = p1[k] - p2[k];
										s += x * x;
									}// k
									double val = s / nCols;
									pDest->add(j, i, val);
								}// j
							}// i
							//
							this->m_rowdist = oDist;
							pRet = pDest;
						}// pDest
					}// ok
				}// pRet
				return (pRet);
			});
		}// getRowDistAsync
		PDistanceMapType get_rows_distances_map(void) {
			PDistanceMapType pRet = this->m_rowdist.get();
			if (pRet == nullptr) {
				pplx::task<PDistanceMapType> tsk = this->getRowDistAsync();
				pRet = tsk.get();
			}
			return (pRet);
		}// get_rows_distances_map
		pplx::task<PDistanceMapType> getColDistAsync(void) {
			return pplx::task<PDistanceMapType>([this]()->PDistanceMapType {
				PDistanceMapType pRet = this->m_coldist.get();
				if (pRet == nullptr) {
					const size_t nRows = this->m_nrows;
					const size_t nCols = this->m_ncols;
					const double *pf = this->m_data.get();
					if ((nRows > 0) && (nCols > 0) && (pf != nullptr)) {
						std::shared_ptr<DistanceMapType> oDist(new DistanceMapType());
						PDistanceMapType pDest = oDist.get();
						if (pDest != nullptr) {
							//
							for (size_t i = 0; i < nCols; ++i) {
								for (size_t j = 0; j < i; ++j) {
									double s = 0;
									for (size_t k = 0; k < nRows; ++k) {
										double x = pf[k * nCols + i] - pf[k * nCols + j];
										s += x * x;
									}// k
									double val = s / nRows;
									pDest->add(j, i, val);
								}// j
							}// i
							 //
							this->m_coldist = oDist;
							pRet = pDest;
						}// pDest
					}// ok
				}// pRet
				return (pRet);
			});
		}// getColistAsync
		PDistanceMapType get_cols_distances_map(void) {
			PDistanceMapType pRet = this->m_coldist.get();
			if (pRet == nullptr) {
				pplx::task<PDistanceMapType> tsk = this->getColDistAsync();
				pRet = tsk.get();
			}
			return (pRet);
		}// get_cols_distances_map
		/////////////////////////////////
		template< typename X>
		static pplx::task<MatDataPtr> createAsync(const size_t nRows, const size_t nCols, const X &pxdata) {
			return pplx::task<MatDataPtr>([nRows, nCols, pxdata]()->MatDataPtr {
				MatDataPtr oRet = std::make_shared<MatData>();
				MatData *pRet = oRet.get();
				if (pRet != nullptr) {
					pRet->initAsync(nRows, nCols, pxdata).then([oRet](bool b) {
						if (b) {
							MatData *px = oRet.get();
							auto t1 = px->getRowDistAsync();
							auto t2 = px->getRowDistAsync();
							auto tr = t1 && t2;
							tr.wait();
						}
						return (oRet);
					});
				}// pRet
				return (oRet);
			});
		}// create
		/////////////////////////////////////
	};// class MatData
	///////////////////////////////////////
}// namespace info
////////////////////////////////////////////
#endif // !__MATDATA_H__

