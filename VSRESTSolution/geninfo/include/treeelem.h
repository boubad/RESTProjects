#pragma once
#ifndef __INDIVTREE_H__
#define __INDIVTREE_H__
/////////////////////////////////////
#include "matdata.h"
#include <map>
//////////////////////////////
namespace info {
////////////////////////////////////
template<typename U,typename DATATYPE, typename DISTANCETYPE>
class TreeItem  {
public:
	using index_type  = U;
	using data_type = DATATYPE;
	using distance_type = DISTANCETYPE;
	//
	using TreeItemType = TreeItem<index_type,data_type,distance_type>;
	using PTreeItemType = TreeItemType *;
	using ints_sizet_map = std::map<index_type,size_t>;
	using ints_vector = std::vector<index_type>;
private:
	size_t m_ncols;
	index_type m_index;
	data_type *m_pcenter;
	TreeItemType *m_pleft;
	TreeItemType *m_pright;
	//
	TreeItem(const TreeItemType &) = delete;
	TreeItemType & operator=(const TreeItemType &) = delete;
public:
	template <typename X>
	TreeItem(index_type aIndex,size_t nCols, const X &pxdata) : 
	m_ncols(nCols),m_index(aIndex),m_pcenter(nullptr),m_pleft(nullptr),m_pright(nullptr) {
	 assert(aIndex != 0); 
	 assert(nCols > 0);
	 data_type *p = new data_type[nCols];
	 assert(p != nullptr);
	 for (size_t i = 0; i < nCols; ++i){
	   p[i] = (data_type)pxdata[i];
	 }
	 this->m_pcenter = p;
	}
	TreeItem(TreeItemType *pLeft, TreeItemType *pRight) :
	m_ncols(0),m_index(0),m_pcenter(nullptr),m_pleft(nullptr),m_pright(nullptr) {
		assert(pLeft != nullptr);
		assert(pRight != nullptr);
		const size_t nCols = pLeft->m_ncols;
		assert(nCols > 0);
		assert(pRight->m_ncols == nCols);
		this->m_ncols = nCols;
		this->m_pleft = pLeft;
		this->m_pright = pRight;
		data_type *p = new data_type[nCols];
		assert(p != nullptr);
		data_type *p1 = pLeft->m_pcenter;
		assert(p1 != nullptr);
		data_type *p2 = pRight->m_pcenter;
		assert(p2 != nullptr);
		for (size_t i = 0; i < nCols; ++i){
		  p[i] = (data_type)((p1[i] + p2[i]) / 2.0);
		}// i
		this->m_pcenter = p;
	} // TreeItem
	virtual ~TreeItem() {
		delete this->m_pcenter;
		delete this->m_pleft;
		delete this->m_pright;
	}
public:
	bool is_leaf(void) const {
		return ((this->m_pleft == nullptr) && (this->m_pright == nullptr));
	}
	const data_type * center(void) const {
		return (this->m_pcenter);
	}
	const TreeItemType *left(void) const {
		return (this->m_pleft);
	}
	const TreeItemType *right(void) const {
		return (this->m_pright);
	}
	void get_map(ints_sizet_map &oMap, const size_t val) const {
		if (this->is_leaf()) {
		  oMap[this->m_index] = val;
		return;
		}
		if (this->m_pleft != nullptr) {
			this->m_pleft->get_map(oMap, val);
		}
		if (this->m_pright != nullptr) {
			this->m_pright->get_map(oMap, val);
		}
	} // get_map
	void get_ids(ints_vector &oVec) const {
		if (this->is_leaf()) {
			oVec.push_back(this->m_index);
			return;
		}
		if (this->m_pleft != nullptr) {
			this->m_pleft->get_ids(oVec);
		}
		if (this->m_pright != nullptr) {
			this->m_pright->get_ids(oVec);
		}
	} // get_map
	distance_type distance(const TreeItemType &other) const {
		const data_type *pc1 = this->m_pcenter;
		assert(pc1 != nullptr);
		const data_type *pc2 = other.m_pcenter;
		assert(pc2 != nullptr);
		distance_type dRes = 0;
		const size_t nCols = this->m_ncols;
		assert(nCols > 0);
		assert(other.m_ncols == nCols);
		for (size_t i = 0; i < nCols; ++i){
		   distance_type x = (distance_type)(pc1[i] - pc2[i]);
		   dRes = (distance_type)(dRes + x * x);
		}// i
		return (dRes);
	} // distance
};
// class TreeItem<U>
/////////////////////////////////////////
}// namespace info
////////////////////////////////////////
#endif // !__INDIVTREE_H__
