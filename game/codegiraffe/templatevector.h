#pragma once

#include "templatearray.h"

/**
 * a custom implementation of a vector, similar to std::vector
 * Ideal for dynamic lists of primitive types and single pointers.
 * written for understandability and convenience. for pure efficiency, use std::vector
 * @WARNING template with virtual types, or types that will be pointed to at your own risk!
 * in those situations, templates of pointers to those types are a much better idea.
 *
 * @author mvaganov@hotmail.com October 2014
 */
template<typename DATA_TYPE>
class TemplateVector : public TemplateArray<DATA_TYPE>
{
protected:
	/** the default size to grow vectors by */
	static const int DEFAULT_ALLOCATION_SIZE = 10;

	/** number of valid elements that the caller thinks we have.. */
	int m_size;

public:
	/** @return the size of the list */
	int size() const { return m_size; }

	/** sets all fields to an initial data state. WARNING: can cause memory leaks if used without care */
	void init() {
		TemplateArray<DATA_TYPE>::init();
		m_size = 0;
	}

	/** cleans up memory */
	void release() {
		TemplateArray<DATA_TYPE>::release();
		m_size = 0;
	}

	/** @return true if the copy finished correctly */
	bool copy(TemplateVector<DATA_TYPE> const & a_vector) {
		if(ensureCapacity(a_vector.m_size)) {
			for(int i = 0; i < a_vector.m_size; ++i) {
				set(i, a_vector.get(i));
			}
			m_size = a_vector.m_size;
			return true;
		}
		return false;
	}

	/** copy constructor */
	TemplateVector(TemplateVector<DATA_TYPE> const & a_vector) {
		init();
		copy(a_vector);
	}

	/** default constructor allocates no list (zero size) */
	TemplateVector(){ init(); }

	/** add an array to this vector */
	void add(int const a_size, DATA_TYPE * const & a_values) {
		ensureCapacity(size()+a_size);
		for(int i = 0; i < a_size; ++i)
			add(a_values[i]);
	}

	/** format constructor */
	TemplateVector(int const a_size, DATA_TYPE const & a_defaultValue) {
		init();
		add(a_size, a_defaultValue);
	}

	/** complete constructor */
	TemplateVector(int const a_size, DATA_TYPE * const & a_defaultValues) {
		init();
		ensureCapacity(a_size);
		for(int i = 0; i < a_size; ++i)
			add(a_defaultValues[i]);
	}

	/** @return the last value in the list */
	DATA_TYPE & getLast() {
		return m_data[m_size-1];
	}

	/** @return the last added value in the list, and lose that value */
	DATA_TYPE & pop() {
		return m_data[--m_size];
	}

	/** @return add an element to the business end of this list structure */
	void push(DATA_TYPE const & a_value) { add(a_value); }

	/**
	 * @param value to add to the list 
	 * @note adding a value may cause memory allocation
	 */
	void add(DATA_TYPE const & a_value) {
		// where am i storing these values?
		// if i don't have a place to store them, i better make one.
		if(!m_data) {
			// make a new list to store numbers in
			allocateToSize(DEFAULT_ALLOCATION_SIZE);
		}
		// if we don't have enough memory allocated for this list
		if(m_size >= m_allocated) {
			// make a bigger list
			allocateToSize(m_allocated*2);
		}
		set(m_size, a_value);
		m_size++;
	}

	/**
	 * @param a_value to add to the list if it isnt in the list already
	 * @return true the index where the element exists
	 * @note see also: insertSorted(), which also supports no duplicates
	 */
	int addNoDuplicates(DATA_TYPE const & a_value) {
		int index = indexOf(a_value);
		if(index < 0) {
			index = m_size;
			add(a_value);
		}
		return index;
	}

	/** @param a_vector a vector to add all the elements from */
	void addVector(TemplateVector<DATA_TYPE> const & a_vector) {
		for(int i = 0; i < a_vector.size(); ++i) {
			add(a_vector.get(i));
		}
	}

	/** 
	 * @param size the user wants the vector to be (chopping off elements)
	 * @return false if could not allocate memory
	 * @note may cause memory allocation if size is bigger than current
	 */
	bool setSize(int const a_size) {
		if(!ensureCapacity(a_size))
			return false;
		m_size = a_size;
		return true;
	}

	/** sets size to zero, but does not deallocate any memory */
	void clear() { setSize(0); }

	/** 
	 * @param a_index is overwritten by the next element, which is 
	 * overwritten by the next element, and so on, till the last element
	 */
	void remove(int const a_index) {
		moveDown(a_index, -1, m_size);
		setSize(m_size-1);
	}

	/** 
	 * @param a_index where to insert a_value. shifts elements in the vector.
	 */
	void insert(int const a_index, DATA_TYPE const & a_value) {
		setSize(m_size+1);
		moveUp(a_index, 1, m_size);
		set(a_index, a_value);
	}

	/** 
	 * @return first element from the list and moves the rest up 
	 * @note removes the first element from the list
	 */
	DATA_TYPE pull() {
		DATA_TYPE value = get(0);
		remove(0);
		return value;
	}

	/** @param a_index is replaced by the last element, then size is reduced. */
	void removeFast(int const a_index) {
		swap(a_index, m_size-1);
		setSize(m_size-1);
	}

	/** @return the index of the first appearance of a_value in this vector. uses == */
	int indexOf(DATA_TYPE const & a_value) const {
		return TemplateArray<DATA_TYPE>::indexOf(a_value, 0, m_size);
	}

	/** @return index of 1st a_value at or after a_startingIndex. uses == */
	int indexOf(DATA_TYPE const & a_value, int const a_startingIndex) const {
		return TemplateArray<DATA_TYPE>::indexOf(a_value, a_startingIndex, m_size);
	}

	/**
	 * will only work correctly if the TemplateVector is sorted.
	 * @return the index of the given value, -1 if the value is not in the list
	 */
	int indexOfWithBinarySearch(DATA_TYPE const & a_value) const {
		if(m_size) {
			return TemplateArray::indexOfWithBinarySearch(a_value, 0, m_size-1);
		}
		return -1;    // failed to find key
	}

	/**
	* uses binary sort to put values in the correct index.
	* @param a_value value to insert in order
	* @return the index where a_value was inserted
	*/
	int insertAsSet(DATA_TYPE const & a_value) { return insertSorted(a_value, false); }

	/**
	 * uses binary sort to put values in the correct index. safe if sorting is always used
	 * @param a_value value to insert in order
	 * @param a_allowDuplicates will not insert duplicates if set to false
	 * @return the index where a_value was inserted
	 */
	int insertSorted(DATA_TYPE const & a_value, bool const a_allowDuplicates) {
		const int a_start = 0, a_size = m_size;
		int index = 0;
		if (a_size == 0 || a_value < m_data[a_start])
			index = a_start;
		else if (m_data[m_size - 1] < a_value) {
			index = m_size;
		} else {
			int imin = a_start, imax = a_start + a_size - 1;
			while (imax >= imin) {
				index = (imin + imax) / 2; // calculate the midpoint for roughly equal partition
				if (m_data[index] == a_value) { // key found at index. 
					if (a_allowDuplicates) { // if there are duplicates, and that is ok, find the end of the duplicates
						while (index + 1 < m_size && m_data[index + 1] == a_value) index++;
						index++; // go one more, so this is appended to the end of the stream of duplicates
					}
					break;
				} else if (m_data[index] < a_value) { // determine which subarray to search
					imin = index + 1; // change min index to search upper subarray
				} else {
					imax = index - 1; // change max index to search lower subarray
				}
			}
			if (imin > imax) {
				if (imin == index+1) index = imin;
				//else if (imax == index - 1) index = imax;
			}
		}
		if(!m_size  || a_allowDuplicates || !(a_value == m_data[index]))
			insert(index, a_value);
		//if (!isSorted()) { printf("bad news... it's not sorted anymore...\n"); }
		return index;

	}

	/**
	 * @param a_value first appearance replaced by last element. breaks if not in list
	 */
	void removeDataFast(DATA_TYPE const & a_value) {
		removeFast(indexOf(a_value));
	}

	/**
	 * @param a_listToExclude removes these elements from *this list
	 * @return true if at least one element was removed
	 */
	bool removeListFast(TemplateVector<DATA_TYPE> const & a_listToExclude) {
		bool aTermWasRemoved = false;
			for(int i = size()-1; i >- 0; --i) {
				for (int e = 0; e < a_listToExclude.size(); ++e) {
					if (a_listToExclude.get(e) == get(i)) {
					removeFast(i);
					aTermWasRemoved = true;
				}
			}
		}
		return aTermWasRemoved;
	}

	/**
	 * @param a_value first appearance is removed. 
	 * @return if data was removed
	 */
	bool removeData(DATA_TYPE const & a_value) {
		int index = indexOf(a_value);
		if(index >= 0) {
			remove(index);
			return true;
		}
		return false;
	}

	/** destructor */
	~TemplateVector() { release(); }

	/** re-defined here because the nature of the size() method has changed */
	void reverse() {
		int half = size() / 2;
		for (int i = 0; i < half; ++i) {
			swap(i, size() - 1 - i);
		}
	}

	bool isSorted() { return TemplateArray<DATA_TYPE>::isSorted(0, size()); }

	void quicksort() { TemplateArray<DATA_TYPE>::quicksort(0, size() - 1); }
};
