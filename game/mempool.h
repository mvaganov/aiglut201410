#pragma once

#include "codegiraffe/templatevectorlist.h"
#include "codegiraffe/templateset.h"

template <typename TYPE>
class MemPool
{
	TemplateSet<TYPE*> freeList;
	TemplateVectorList<TYPE> list;
public:
	MemPool(){}

	TYPE * add() {
		if (freeList.size() == 0) {
			list.add();
			return &list.getLast();
		}
		return freeList.pop();
	}
	TYPE & add(TYPE const & a_data) {
		TYPE * data = add();
		*data = a_data;
		return *data;
	}
	void markFree(TYPE* data) {
		int index = list.indexOf(data);
		// being asked to mark an element that is not actually in here is a problem.
		if (index == -1) {
			int i = 0; i = 1 / i;
		}
		if (index == list.size() - 1) {
			list.removeFast(list.size() - 1);
		} else {
			freeList.add(data);
			if (countUsed() == 0) {
				clear();
			}
		}
	}
	void clear() {
		list.clear();
		freeList.clear();
	}
	int countUsed() const { return list.size() - freeList.size(); }
	int size() const { return list.size(); }
	int countFree() const { return freeList.size(); }
	TYPE& get(const int index) { return list[index]; }
	TYPE& operator[](const int i) { return list[i]; }
	const TYPE& operator[](const int i) const { return list[i]; }
	bool isMarkedFree(TYPE * data) const { return freeList.has(data); }
};