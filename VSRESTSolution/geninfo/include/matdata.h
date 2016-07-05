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
#include "distance_map.h"
//////////////////////
namespace info {
	////////////////////////////////////////
	template <typename INDEXTYPE, typename FLOATTYPE, typename DISTANCETYPE, typename STRINGTYPE>
	class MatData {
	public:
		using index_type = INDEXTYPE;
		using data_type = FLOATTYPE;
		using distance_type = DISTANCETYPE;
		using string_type = STRINGTYPE;
		//
		using strings_vector = std::vector<string_type>;
		using data_vector = std::vector<data_type>;
		using distancemap_type = DistanceMap<index_type, distance_type>;
		using matdata_type = MatData<index_type, data_type, distance_type, string_type>;
		using matdata_type_ptr = std::shared_ptr<matdata_type>;
	private:
		size_t m_nrows;
		size_t m_ncols;
		std::shared_ptr<data_type> m_data;
		strings_vector m_rownames;
		strings_vector m_colnames;
		std::shared_ptr<distancemap_type> m_rowdist;
		std::shared_ptr<distancemap_type> m_coldist;
	public:
		MatData() : m_nrows(0), m_ncols(0) {
		}
		MatData(const matdata_type &other) :
			m_nrows(other.m_nrows), m_ncols(other.m_ncols), m_data(
				other.m_data), m_rownames(other.m_rownames), m_colnames(
					other.m_colnames), m_rowdist(other.m_rowdist), m_coldist(
						other.m_coldist) {
		}
		matdata_type & operator=(const matdata_type &other) {
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
			return ((this->m_nrows > 0) && (this->m_ncols > 0) && (this->m_data.get() != nullptr));
		}	// is_valid
		size_t rows(void) const {
			return (this->m_nrows);
		}	// rows
		size_t cols(void) const {
			return (this->m_ncols);
		}	// cols
		data_type get_data(const size_t irow, const size_t icol) const {
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
		string_type get_row_name(const size_t irow) const {
			const strings_vector &v = this->m_rownames;
			assert((irow < this->m_nrows) && (irow < v.size()));
			return (v[irow]);
		}	// row_name
		string_type get_col_name(const size_t icol) const {
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
			if ((pColNames != nullptr) && (pColNames->size() >= nCols)) {
				const strings_vector &v = *pColNames;
				for (size_t i = 0; i < nCols; ++i) {
					string_type s = v[i];
					colnames[i] = s;
				}	// i
			}
			this->m_rownames = rownames;
			this->m_colnames = colnames;
			return (true);
		}	// setnames
	public:
		template<typename X>
		bool init(size_t nRows, size_t nCols, const X &pxdata,
			const strings_vector *pRowNames = nullptr,
			const strings_vector *pColNames = nullptr) {
			bool bRet = false;
			if ((nCols < 1) || (nRows < 1)) {
				return (bRet);
			}
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
						double x = (double)pxdata[j*nCols + i];
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
						pf[j * nCols + i] = (data_type)((pxdata[j*nCols + i] - vmin) / delta);
					}	// j
				}	// i
				this->m_nrows = nRows;
				this->m_ncols = nCols;
				this->m_data = ff;
				this->m_rowdist.reset();
				this->m_coldist.reset();
				this->set_names(pRowNames, pColNames);
				bRet = true;
			}
			catch (...) {}
			return (bRet);
		}	// init
		distancemap_type *get_rows_distances_map(void) {
			distancemap_type  *pRet = this->m_rowdist.get();
			if (pRet == nullptr) {
				this->getRowDist();
				pRet = this->m_rowdist.get();
			}
			return (pRet);
		}// get_rows_distances_map
		distancemap_type  *get_cols_distances_map(void) {
			distancemap_type  *pRet = this->m_coldist.get();
			if (pRet == nullptr) {
				this->getColDist();
				pRet = this->m_coldist.get();
			}
			return (pRet);
		}// get_cols_distances_map
		template<typename X>
		static matdata_type_ptr create(size_t nRows, size_t nCols, const X &pxdata,
			const strings_vector *pRowNames = nullptr,
			const strings_vector *pColNames = nullptr) {
			matdata_type_ptr oRet = std::make_shared<matdata_type>();
			matdata_type *pRet = oRet.get();
			if (pRet != nullptr) {
				if (pRet->init(nRows, nCols, pxdata, pRowNames, pColNames)) {
					pRet->getRowDist();
					pRet->getColDist();
				}
			}// pRet
			return (oRet);
		}		// create
		template<typename X>
		static std::future<matdata_type_ptr> createAsync(size_t nRows, size_t nCols, const X &pxdata,
			const strings_vector *pRowNames = nullptr,
			const strings_vector *pColNames = nullptr) {
			return std::async([nRows, nCols, pxdata, pRowNames, pColNames]()->matdata_type_ptr {
				return MatDataType::create(nRows, nCols, pxdata, pRowNames, pColNames);
			});
		}		// createAsync
				/////////////////////////////////////
	protected:
		void getRowDist(void) {
			const size_t nRows = this->m_nrows;
			const size_t nCols = this->m_ncols;
			const data_type *pf = this->m_data.get();
			if ((nRows > 0) && (nCols > 0) && (pf != nullptr)) {
				DistanceMap<size_t, float> oDist;
				for (size_t i = 0; i < nRows; ++i) {
					const data_type *p1 = pf + i * nCols;
					for (size_t j = 0; j < i; ++j) {
						const data_type *p2 = pf + (j * nCols);
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
				this->m_rowdist = std::make_shared<distancemap_type>();
				distancemap_type *pRet = this->m_rowdist.get();
				if (pRet != nullptr) {
					oDist.recode(*pRet);
				}
			}
		}	// getRowDistAsync
		void getColDist(void) {
			const size_t nRows = this->m_nrows;
			const size_t nCols = this->m_ncols;
			const data_type *pf = this->m_data.get();
			if ((nRows > 0) && (nCols > 0) && (pf != nullptr)) {
				DistanceMap<size_t, float> oDist;
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
				this->m_coldist = std::make_shared<distancemap_type>();
				distancemap_type *pRet = this->m_coldist.get();
				if (pRet != nullptr) {
					oDist.recode(*pRet);
				}
			}	 //
		}							 // getColistAsync
	};
	// class MatData
	///////////////////////////////////////
}// namespace info
 ////////////////////////////////////////////
#endif // !__MATDATA_H__
