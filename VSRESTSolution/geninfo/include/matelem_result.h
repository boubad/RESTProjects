#pragma once
#ifndef __MATELEM_RESULT_H__
#define __MATELEM_RESULT_H__
//////////////////////////
#include "gendefs.h"
///////////////////////////////
namespace info {
	//////////////////////////////////
	template <typename INDEXTYPE, typename CRITERIATYPE>
	class MatElemResult {
	public:
		using index_type = INDEXTYPE;
		using criteria_type = CRITERIATYPE;
		using matelemresult_type = MatElemResult<index_type, criteria_type>;
		using matelemresult_type_ptr = std::shared_ptr<matelemresult_type>;
	private:
		DispositionType m_disp;
		criteria_type m_crit;
		size_t m_nsize;
		std::unique_ptr<index_type> m_index;
	public:
		MatElemResult(DispositionType disp = DispositionType::dispUnknown) :m_disp(disp), m_crit(0),m_nsize(0) {}
		template <typename U, typename X>
		MatElemResult(const DispositionType disp, const X crit, const size_t n,const X *pInds) :
			m_disp(disp), m_crit((criteria_type)crit),m_nsize(n),m_index(new index_type[n]) {
			index_type *p = this->m_index.get();
			if ((p != nullptr) && (pInds != nullptr)) {
				for (size_t i = 0; i < n; ++i) {
					p[i] = (index_type)pInds[i];
				}// i
			}// p
		}// MatElemRezsult
		template <typename U, typename X>
		MatElemResult(const DispositionType disp, const X crit, const std::vector<U> &oInds) :
			m_disp(disp), m_crit((criteria_type)crit), m_nsize(oInds.size()), m_index(new index_type[oInds.size()]) {
			index_type *p = this->m_index.get();
			if (p != nullptr) {
				for (size_t i = 0; i < n; ++i) {
					p[i] = (index_type)oInds[i];
				}// i
			}// p
		}// MatElemRezsult
		MatElemResult(const matelemresult_type  &other):m_disp(other.m_disp),m_crit(other.m_crit),
		m_nsize(0){
			const size_t n = other.m_nsize;
			const index_type *pSrc = other.m_index.get();
			if ((pSrc != nullptr) && (n > 0)) {
				this->m_index.reset(new index_type[n]);
				index_type *pDest = this->m_index.get();
				if (pDest != nullptr) {
					this->m_nsize = n;
					for (size_t i = 0; i < n; ++i) {
						pDest[i] = pSrc[i];
					}
				}// pDest
			}
		}
		matelemresult_type  & operator=(const matelemresult_type  &other) {
			if (this != &other) {
				this->m_nsize = 0;
				this->m_index.reset();
				this->m_disp = other.m_disp;
				this->m_crit = other.m_crit;
				const size_t n = other.m_nsize;
				const index_type *pSrc = other.m_index.get();
				if ((pSrc != nullptr) && (n > 0)) {
					this->m_index.reset(new index_type[n]);
					index_type *pDest = this->m_index.get();
					if (pDest != nullptr) {
						this->m_nsize = n;
						for (size_t i = 0; i < n; ++i) {
							pDest[i] = pSrc[i];
						}
					}// pDest
				}
			}
			return (*this);
		}
		virtual ~MatElemResult() {}
		bool operator<(const matelemresult_type &other) const {
			return (this->m_crit < other.m_crit);
		}// operator<
		bool operator==(const matelemresult_type &other) const {
			return (this->m_crit == other.m_crit);
		}// operator==
		bool  add(const matelemresult_type &other) {
			bool bRet = false;
			if ((this->m_crit > 0) && (other.m_crit > 0)) {
				if (other.m_crit < this->m_crit) {
					(*this) = other;
					bRet = true;
				}
			}
			else if ((this->m_crit == 0) && (other.m_crit > 0)) {
				(*this) = other;
				bRet = true;
			}
			return (bRet);
		}// add
		matelemresult_type & operator+=(const matelemresult_type &other) {
			this->add(other);
			return (*this);
		}
		std::ostream & write_to(std::ostream &os) const {
			if (this->m_disp == DispositionType::dispCol) {
				os << "COLS ";
			}
			else if (this->m_disp == DispositionType::dispRow) {
				os << "ROWS ";
			}
			os << "CRIT: " << this->m_crit << " [";
			const size_t n = this->m_nsize;
			const index_type *p = this->m_index.get();
			if ((n > 0) && (p != nullptr)) {
				for (size_t i = 0; i < n; ++i) {
					if (i != 0) {
						os << ", ";
					}
					os << p[i];
				}// i
			}
			os << "]";
			return (os);
		}// write_to
		std::wostream & write_to(std::wostream &os) const {
			if (this->m_disp == DispositionType::dispCol) {
				os << L"COLS ";
			}
			else if (this->m_disp == DispositionType::dispRow) {
				os << L"ROWS ";
			}
			os << L"CRIT: " << this->m_crit << L" [";
			const size_t n = this->m_nsize;
			const index_type *p = this->m_index.get();
			if ((n > 0) && (p != nullptr)) {
				for (size_t i = 0; i < n; ++i) {
					if (i != 0) {
						os << L", ";
					}
					os << p[i];
				}// i
			}
			os << L"]";
			return (os);
		}// write_to
	public:
		bool is_valid(void) const {
			return ((this->m_nsize > 0) && (this->m_index.get() != nullptr));
		}// is_valid
		size_t size(void) const {
			return (this->m_nsize);
		}
		DispositionType disposition(void) const {
			return (this->m_disp);
		}
		void disposition(DispositionType disp) {
			this->m_disp = disp;
		}
		criteria_type criteria(void) const {
			return (this->m_crit);
		}
		template <typename X>
		void criteria(const X c) {
			this->m_crit = (criteria_type)c;
		}
		template <typename X>
		void indexes(std::vector<X> &inds) const {
			const size_t n = this->m_nsize;
			inds.resize(n);
			const index_type *p = this->m_index.get();
			if ((n > 0) && (p != nullptr)) {
				for (size_t i = 0; i < n; ++i) {
					inds[i] = (X)p[i];
				}
			}
		}
		template <typename X>
		void indexes( X *pInds) const {
			assert(pInds != nullptr);
			const size_t n = this->m_nsize;
			const index_type *p = this->m_index.get();
			if ((n > 0) && (p != nullptr)) {
				for (size_t i = 0; i < n; ++i) {
					pInds[i] = (X)p[i];
				}
			}
		}
		template <typename U>
		bool permute(const U i, const U j) {
			const size_t n = this->m_nsize;
			index_type *p = this->m_index.get();
			size_t ii = (size_t)i;
			size_t jj = (size_t)j;
			if ((p == nullptr) || (ii >= n) || (jj >= n)) {
				return (false);
			}
			index_type t = p[ii];
			p[ii] = p[jj];
			p[jj] = t;
			return (true);
		}// permute
		template <typename U, typename X>
		void set(const  X crit, const std::vector<U> &oInds) {
			this->m_crit = (criteria_type)crit;
			const size_t n = oInds.size();
			this->m_nsize = n;
			this->m_index.reset(new index_type[n]);
			index_type *p = this->m_index.get();
			if (p != nullptr) {
				for (size_t i = 0; i < n; ++i) {
					p[i] = (index_type)oInds[i];
				}// i
			}// p
		}// set
		template <typename U, typename X>
		void set(const  X crit, const size_t n, const U *pInds) {
			this->m_crit = (criteria_type)crit;
			const size_t n = oInds.size();
			this->m_nsize = n;
			this->m_index.reset(new index_type[n]);
			index_type *p = this->m_index.get();
			if ((p != nullptr) && (pInds != nullptr)) {
				for (size_t i = 0; i < n; ++i) {
					p[i] = (index_type)pInds[i];
				}// i
			}// p
		}// set
		template <typename U, typename X>
		static matelemresult_type_ptr create(const DispositionType disp, const X crit, const std::vector<U> &oInds) {
			return (std::make_shared<matelemresult_type>(disp, crit, oInds));
		}// create
		template <typename U, typename X>
		static matelemresult_type_ptr create(const DispositionType disp, const X crit, const size_t n,const X *pInds) {
			return (std::make_shared<matelemresult_type>(disp, crit, n ,pInds));
		}// create
	};
	//////////////////////////////
	template <typename INDEXTYPE, typename CRITERIATYPE>
	class MatElemFunction {
	public:
		using index_type = INDEXTYPE;
		using criteria_type = CRITERIATYPE;
		using matelemresult_type = MatElemResult<index_type, criteria_type>;
		using matelemresult_type_ptr = std::shared_ptr<matelemresult_type>;
		using matelemfunction_type = MatElemFunction<index_type, criteria_type>;
	public:
		MatElemFunction() {}
		MatElemFunction(const matelemfunction_type &) {}
		matelemfunction_type & operator=(const matelemfunction_type &) {
			return (*this);
		}
		virtual ~MatElemFunction() {}
	public:
		virtual void perform(const matelemresult_type_ptr r) const {
			// do nothing here
		}
		void operator()(const matelemresult_type_ptr r) const {
			this->perform(r);
		}// operator()
	}; // class MatElemFunction
	   ////////////////////////
} // namespace info
////////////////////////////////////////////
template <typename INDEXTYPE, typename CRITERIATYPE>
std::ostream & operator<<(std::ostream &os, const info::MatElemResult<INDEXTYPE, CRITERIATYPE> &p) {
	return (p.write_to(os));
}// operator<<
template <typename INDEXTYPE,typename CRITERIATYPE>
info::MatElemResult<INDEXTYPE, CRITERIATYPE> operator+(const info::MatElemResult<INDEXTYPE, CRITERIATYPE> &p1,
	const info::MatElemResult<INDEXTYPE, CRITERIATYPE> &p2) {
	info::MatElemResult<INDEXTYPE, CRITERIATYPE> temp(p1);
	(void)temp.add(p2);
	return (temp);
}// operator+
  /////////////////////////////////
#endif  //__MATELEM_RESULT_H__
