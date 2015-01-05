#pragma once
#include "templatevector.h"

template<typename TYPE>
class TemplateSet : protected TemplateVector<TYPE> {
public:
	TemplateSet() :TemplateVector<TYPE>(){}

	int size() const { return TemplateVector<TYPE>::size(); }
	bool add(TYPE const & a_value) { return TemplateVector<TYPE>::insertSorted(a_value, false) != -1; }
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
	TYPE * getRawList() { return TemplateVector<TYPE>::getRawList(); }
	const TYPE * getRawListConst() const { return TemplateVector<TYPE>::getRawListConst(); }
	TYPE & getLast() { return TemplateVector<TYPE>::getLast(); }
	bool operator==(TemplateSet<TYPE> const & set) const { return TemplateVector<TYPE>::operator==(set); }
	TemplateSet & operator=(TemplateSet<TYPE> const & other) { TemplateVector<TYPE>::operator=(other); return *this; }
	bool containsAll(const TYPE * subset, const int subsetCount) const {
		bool found = false;
		for (int i = 0; i < subsetCount; ++i) {
			found = indexOf(subset[i]) >= 0;
			if (!found) break;
		}
		return found;
	}
	void addVector(TemplateVector<TYPE> & a_vector) {
		for (int i = 0; i < a_vector.size(); ++i) {
			add(a_vector[i]);
		}
	}
};