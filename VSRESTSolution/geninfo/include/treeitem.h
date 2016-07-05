#pragma once
#ifndef __TREEELEM_H__
#define __TREEELEM_H___
/////////////////////////////////////
#include "indiv.h"
///////////////////////////////
namespace info {
	////////////////////////////////////
	template<typename INDEXTYPE, typename DATATYPE, typename DISTANCETYPE, typename STRINGTYPE>
	class TreeItem : public Indiv<INDEXTYPE, DATATYPE, STRINGTYPE> {
	public:
		using index_type = INDEXTYPE;
		using data_type = DATATYPE;
		using distance_type = DISTANCETYPE;
		using string_type = STRINGTYPE;
		//
		using indiv_type = Indiv<index_type, data_type, string_type>;
		using cancellableflag_type = CancellableObject::cancellableflag_type;
		using treeitem_type = TreeItem<index_type, data_type, distance_type, string_type>;
		using treeitem_type_ptr = std::shared_ptr<treeitem_type>;
	private:
		treeitem_type_ptr m_left;
		treeitem_type_ptr m_right;
	public:
		TreeItem(const treeitem_type &other) :indiv_type(other), m_left(other.m_left),
			m_right(other.m_right) {}
		treeitem_type & operator=(const treeitem_type &other) {
			if (this != &other) {
				index_type::operator=(other);
				this->m_left = other.m_left;
				this->m_right = other.m_right;
			}
			return (*this);
		}
		template <typename X>
		TreeItem(index_type aIndex, size_t nCols, const X &pxdata,
			const string_type &sSigle = string_type(),
			cancellableflag_type *pf = nullptr, Backgrounder *pb = nullptr) :
			indiv_type(aIndex, nCols, pxdata, sSigle, pf, pb) {
			assert(aIndex != 0);
		}
		TreeItem(treeitem_type_ptr left, treeitem_type_ptr right,
			cancellableflag_type *pf = nullptr, Backgrounder *pb = nullptr) :indiv_type(pf, pb), m_left(left), m_right(right) {
			treeitem_type *pLeft = left.get();
			assert(pLeft != nullptr);
			treeitem_type *pRight = right.get();
			assert(pRight != nullptr);
			assert(pLeft != pRight);
			const size_t nCols = pLeft->size();
			assert(pRight->size() == nCols);
			this->m_ncols = nCols;
			this->m_center.reset(new data_type[nCols]);
			data_type *p = this->m_center.get();
			assert(p != nullptr);
			const data_type *p1 = pLeft->center();
			const data_type *p2 = pRight->center();
			for (size_t i = 0; i < nCols; ++i) {
				p[i] = (data_type)((p1[i] + p2[i]) / 2.0);
			}// i
		} // TreeItem
		virtual ~TreeItem() {
		}
	public:
		bool is_leaf(void) const {
			return ((this->m_left.get() == nullptr) || (this->m_right.get() == nullptr));
		}
		treeitem_type *left(void) const {
			return (this->m_left.get());
		}
		treeitem_type *right(void) const {
			return (this->m_right.get());
		}
		template <typename XK, typename XV>
		void get_map(std::map<XK, XV> &oMap, const XV val) {
			if (this->is_cancellation_requested()) {
				return;
			}
			if (this->is_leaf()) {
				XK key = (XK)this->id();
				oMap[key] = val;
				return;
			}
			this->left()->get_map(oMap, val);
			this->right()->get_map(oMap, val);
		} // get_map
		template <typename X>
		void get_ids(std::vector<X> &oVec) {
			if (this->is_cancellation_requested()) {
				return;
			}
			if (this->is_leaf()) {
				X ind = (X)this->id();
				oVec.push_back(ind);
				return;
			}
			this->left()->get_ids(oVec);
			this->right()->get_ids(oVec);
		} // get_map
	};
	// class TreeItem<U>
	/////////////////////////////////////////
}// namespace info
 ////////////////////////////////////////
#endif // !__INDIVTREE_H__
