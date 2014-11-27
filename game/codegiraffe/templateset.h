#pragma once
#include "templatevector.h"

template<typename TYPE>
class TemplateSet : protected TemplateVector<TYPE> {
public:
	TemplateSet() :TemplateVector<TYPE>(){}

	int size() const { return TemplateVector<TYPE>::size(); }
	void add(TYPE const & a_value) { TemplateVector<TYPE>::insertSorted(a_value, false); }
	int indexOf(TYPE const & a_value) const {
		return TemplateVector<TYPE>::indexOfWithBinarySearch(a_value);
	}
	bool removeData(TYPE const & a_value) {
		int index = indexOf(a_value);
		if (index >= 0) {
			remove(index);
			return true;
		}
		return false;
	}
	void remove(const int index) { TemplateVector<TYPE>::remove(index); }
	bool has(TYPE const & a_value) const { return indexOf(a_value) >= 0; }
	TYPE & operator[](const int index) { return TemplateVector<TYPE>::operator[](index); }
	const TYPE & operator[](const int index) const { return TemplateVector<TYPE>::operator[](index); }
	void clear() { TemplateVector<TYPE>::clear(); }
	void release() { TemplateVector<TYPE>::release(); }
	bool ensureCapacity(const int a_size) { return TemplateVector<TYPE>::ensureCapacity(a_size); }
	TemplateVector<TYPE> * asVector() { return (TemplateVector<TYPE>*)this; }
	TYPE & pop() { return TemplateVector<TYPE>::pop(); }
};