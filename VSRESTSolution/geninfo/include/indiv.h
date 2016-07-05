#pragma once
#ifndef __INDIV_H__
#define __INDIV_H__
////////////////////////////////////
#include "cancellable_object.h"
////////////////////////////
namespace info {
	///////////////////////////////
	template <typename INDEXTYPE, typename DATATYPE, typename STRINGTYPE>
	class Indiv : public CancellableObject {
	public:
		using string_type = STRINGTYPE;
		using index_type = INDEXTYPE;
		using data_type = DATATYPE;
		using cancellableflag_type = CancellableObject::cancellableflag_type;
		//
		using indiv_type = Indiv<index_type, data_type, string_type>;
	protected:
		index_type m_index;
		size_t m_ncols;
		string_type m_sigle;
		std::shared_ptr<data_type> m_center;
	public:
		Indiv(cancellableflag_type *pf = nullptr,Backgrounder *pb = nullptr) :CancellableObject(pf,pb), m_index(0), m_ncols(0) {}
		Indiv(const indiv_type &other) :CancellableObject(other), m_index(other.m_index), m_ncols(other.m_ncols),
			m_sigle(other.m_sigle), m_center(other.m_center) {}
		index_type & operator=(const index_type &other) {
			if (this != &other) {
				CancellableObject::operator=(other);
				this->m_index = other.m_index;
				this->m_ncols = other.m_ncols;
				this->m_sigle = other.m_sigle;
				this->m_center = other.m_center;
			}
			return (*this);
		}
		template <typename X>
		Indiv(index_type aIndex, size_t nCols, const X &pxdata,
			const string_type &sSigle = string_type(),
			cancellableflag_type *pf = nullptr,Backgrounder *pb = nullptr) :CancellableObject(pf,pb),
			m_index(aIndex), m_ncols(nCols), m_sigle(sSigle), m_center(new data_type[nCols]) {
			assert(nCols > 0);
			data_type *p = this->m_center.get();
			assert(p != nullptr);
			for (size_t i = 0; i < nCols; ++i) {
				p[i] = (data_type)pxdata[i];
			}// i
		}// Indiv
		virtual ~Indiv() {
		}
		bool is_valid(void) const {
			return ((this->m_ncols > 0) && (this->m_center.get() != nullptr));
		}// is_valid
		index_type id(void) const {
			return (this->m_index);
		}
		size_t size(void) const {
			assert(this->m_ncols > 0);
			return (this->m_ncols);
		}
		string_type sigle(void) const {
			return (this->m_sigle);
		}
		const data_type *center(void) const {
			assert(this->m_center.get() != nullptr);
			return (this->m_center.get());
		}
		template <typename XU, typename XA, typename XS, typename F>
		void distance(const Indiv<XU, XA, XS> &other, F &res) {
			if (this->is_cancellation_requested()) {
				return;
			}
			if ((!this->is_valid()) || (!other.is_valid())) {
				return;
			}
			const size_t n = this->size();
			assert(other.size() >= n);
			const data_type *p1 = this->center();
			const XA *p2 = other.center();
			double d = 0;
			for (size_t i = 0; i < n; ++i) {
				double x = (double)p1[i] - (double)p2[i];
				d += x * x;
			}// i
			res = (F)d;
		}// distance
	};// class Indiv<INDEXTYPE,DATATYPE,STRINGTYPE>
	  //////////////////////////////////
}// namespace info
 ////////////////////////////////////
#endif // __INDIV_H__
