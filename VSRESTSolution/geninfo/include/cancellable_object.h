#pragma once
#ifndef __CANCELLABLE_OBJECT_H__
#define __CANCELLABLE_OBJECT_H__
///////////////////////////////
#include "active_object.h"
////////////////////////////
namespace info {
	//////////////////////
	class CancellableObject {
	public:
		using cancellableflag_type = std::atomic<bool>;
		using Operation = Backgrounder::Operation;
	private:
		cancellableflag_type *m_pcancel;
		Backgrounder *m_pback;
	public:
		CancellableObject(cancellableflag_type *pf = nullptr,
			Backgrounder *pBack = nullptr) :m_pcancel(pf) {}
		CancellableObject(const CancellableObject &other) :m_pcancel(other.m_pcancel),
			m_pback(other.m_pback) {}
		CancellableObject & operator=(const CancellableObject &other) {
			if (this != &other) {
				this->m_pcancel = other.m_pcancel;
				this->m_pback = other.m_pback;
			}
			return (*this);
		}
		virtual ~CancellableObject() {}
		bool is_cancellation_requested(void) {
			return ((this->m_pcancel != nullptr) && this->m_pcancel->load());
		}//is_cancellation_requested
		cancellableflag_type *get_cancellable_flag(void) const {
			return (this->m_pcancel);
		}// get_cancellable_flag
		void set_cancellable_flag(cancellableflag_type *pf) {
			this->m_pcancel = pf;
		}
		bool has_backgrounder(void) {
			return (this->m_pback != nullptr);
		}
		Backgrounder *get_backgrounder(void) {
			return (this->m_pback);
		}
		void set_backgrounder(Backgrounder *pb) {
			this->m_pback = pb;
		}
		void send(Operation op) {
			if (this->m_pback != nullptr) {
				this->m_pback->send(op);
			}
		}// send
	}; // class CancellableObject
	   ////////////////////////
}// namespace info
 ////////////////////////////////
#endif  // __CANCELLABLE_OBJECT_H__
