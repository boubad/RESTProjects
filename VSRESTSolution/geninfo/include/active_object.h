/*
* activeobject.h
*
*  Created on: 31 mai 2016
*      Author: boubad
*/

#ifndef ACTIVEOBJECT_H_
#define ACTIVEOBJECT_H_
//////////////////////////
#include "gendefs.h"
/////////////////////////////
namespace info {
	////////////////////////
	class Backgrounder;
	////////////////////////
	class Active {
		friend class Backgrounder;
	public:
		using Operation = std::function<void()>;
	private:
		std::atomic<bool> done;
		std::mutex qlock;
		std::queue<Operation> ops_queue;
		std::condition_variable empty;
		std::thread runnable;
	private:
		Active(const Active&) = delete;
		Active& operator=(const Active&) = delete;
	private:
		Active() : done(false) {
			this->runnable = std::thread([this]() {
				while (!this->done.load()) {
					try {
						std::unique_lock<std::mutex> lck(this->qlock);
						this->empty.wait(lck, [this] {return !this->ops_queue.empty(); });
						Operation op = this->ops_queue.front();
						this->ops_queue.pop();
						op();
					}
					catch (...) {}
				} // while
			});
		} // Active
	public:
		~Active() {
			// Schedule a No-Op runnable to flush the dispatch queue
			this->send([this]() {this->done.store(true); });
			runnable.join();
		} // run
		void send(Operation msg_) {
			std::lock_guard<std::mutex> guard(qlock);
			ops_queue.push(msg_);
			empty.notify_all();
		} // send
	};
	// class ActiveObject
	//////////////////////////////////
	class Backgrounder {
	public:
		using Operation = Active::Operation;
	private:
		std::unique_ptr<Active> _active;
	private:
		Backgrounder(const Backgrounder &) = delete;
		Backgrounder & operator=(const Backgrounder &) = delete;
	public:
		Backgrounder() {
			this->_active.reset(new Active());
			assert(this->_active.get() != nullptr);
		}
		virtual ~Backgrounder() {
		}
		void send(Operation op) {
			this->_active->send(op);
		}	 // send
	};
	// class Backgrounder
	//////////////////////////////
}// namespace info
 ///////////////////////
#endif /* ACTIVEOBJECT_H_ */
