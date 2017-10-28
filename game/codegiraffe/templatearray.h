#pragma once

/**
 * a custom implementation of an array, similar to std::array
 *
 * written for understandability and convenience. for pure efficiency, use std::array
 *
 * @author mvaganov@hotmail.com November 2014
 */
template<typename DATA_TYPE>
class TemplateArray {
protected:
	/** pointer to the allocated data */
	DATA_TYPE * m_data;

	/** actual number of allocated elements that we can use */
	int m_allocated;
public:
	/** sets all fields to an initial data state. WARNING: can cause memory leaks if used without care */
	void init() {
		m_data = nullptr;
		m_allocated = 0;
	}

	/** cause a crash if the index is out of bounds. */
	void validateIndex(const int a_index) const { if (a_index < 0 || a_index >= size()) { int i = 0; i = 1 / i; } }

	/** @return the size of the list */
	int const size() const { return m_allocated; }
	
	/** @return value from the list at given index */
	DATA_TYPE get(int const a_index) const {
		validateIndex(a_index < 0); return m_data[a_index]; // (*this)[a_index] should also work as long as m_data is the first property
	}

	/** @return value from the list at given index explicitly by reference */
	DATA_TYPE & getByRef(int const a_index) {
		validateIndex(a_index); return m_data[a_index];
	}

	/** use for const TemplateArray objects */
	const DATA_TYPE & getByRef(int const a_index) const {
		validateIndex(a_index); return m_data[a_index];
	}

	/** simple mutator sets a value in the list */
	void set(int const a_index, DATA_TYPE const a_value) {
		// complex types must overload DATA_TYPE & operator=(const DATA_TYPE &)
		m_data[a_index] = a_value;
	}

	/** cleans up memory */
	void release() {
		if(m_data) {
			delete [] m_data;
			init();
		}
	}

	/** proxy for getByRef */
	DATA_TYPE & operator[](int const a_index) { return getByRef(a_index); }

	/** use for const TemplateArray objects */
	const DATA_TYPE & operator[](int const a_index) const { return getByRef(a_index); }

	/** deep copy with assignment operator, so allocated have their own references */
	TemplateArray<DATA_TYPE> & operator=(TemplateArray<DATA_TYPE> const & src) { copy(src); return *this; }

	/** @return true if this array and the given array contain the same data in the same order */
	bool operator==(TemplateArray<DATA_TYPE> const & arr) const {
		if (size() != arr.size()) return false;
		for (int i = 0; i < size(); ++i) {
			if (get(i) != arr.get(i)) return false;
		}
		return true;
	}

	/** used for deleting elements */
	void moveDown(int const a_from, int const a_offset, int a_last) {
		for(int i = a_from-a_offset; i < a_last; ++i)
			set(i+a_offset, get(i));
	}
	/** used for adding elements */
	void moveUp(int const a_from, int const a_offset, int a_last) {
		for(int i = a_last-a_offset-1; i >= a_from; --i)
			set(i+a_offset, get(i));
	}

public:

	/** @param a_size reallocate the vector to this size */
	bool allocateToSize(int const a_size) {
		// reallocate a new list with the given size
		DATA_TYPE * newList = new DATA_TYPE[a_size];
		// if the list could not allocate, fail...
		if(!newList) { return false; }
		// the temp list is the one we will keep, while the old list will be dropped.
		DATA_TYPE * oldList = m_data;
		// swap done here so set(index, value) can be called instead of the equals operator
		m_data = newList;
		// if there is old data
		if(oldList) {
			// when copying old data, make sure no over-writes happen.
			int smallestSize = m_allocated<a_size?m_allocated:a_size;
			// fill the new list with the old data
			for(int i = 0; i < smallestSize; ++i) {
				set(i, oldList[i]);
			}
			// get rid of the old list (so we can maybe use the memory later)
			delete [] oldList;
		}
		// mark the new allocated size (held size of oldList)
		m_allocated = a_size;
		return true;
	}

	/** note: this method is memory intesive, and should not be in any inner loops... */
	void add(DATA_TYPE const & a_value) {
		allocateToSize(size()+1);
		set(size()-1, a_value);
	}

	/** note: this method is memory intesive, and should not be in any inner loops... */
	DATA_TYPE * add() {
		allocateToSize(size()+1);
		return &get(size()-1);
	}

	/** @return true if vector allocated this size */
	bool ensureCapacity(int const a_size) {
		if(a_size && m_allocated < a_size) {
			return allocateToSize(a_size);
		}
		return true;
	}

	/** @return true of the copy finished correctly */
	bool copy(TemplateArray<DATA_TYPE> const & a_array) {
		if(m_allocated != a_array.m_allocated) {
			release();
			allocateToSize(a_array.m_allocated);
		}
		for(int i = 0; i < a_array.m_allocated; ++i) {
			set(i, a_array.get(i));
		}
		return false;
	}

	/** copy constructor */
	TemplateArray(TemplateArray<DATA_TYPE> const & a_array) {
		init();
		copy(a_array);
	}

	/** default constructor allocates no list (zero size) */
	TemplateArray(){ init(); }


	/** size constructor */
	TemplateArray(int const a_size) {
		init();
		ensureCapacity(a_size);
	}

	/** format constructor */
	TemplateArray(int const a_size, const DATA_TYPE & a_defaultValue) {
		init();
		ensureCapacity(a_size);
		for(int i = 0; i < a_size; ++i)
			set(i, a_defaultValue);
	}

	/** complete constructor */
	TemplateArray(int const a_size, const DATA_TYPE * const & a_defaultValues) {
		init();
		ensureCapacity(a_size);
		for (int i = 0; i < a_size; ++i)
			set(i, a_defaultValues[i]);
	}

	~TemplateArray(){ release(); }
	
	/** @return the last value in the list */
	DATA_TYPE & getLast() { return get(size()-1); }

	/**
	 * @return the raw pointer to the data... 
	 * @note dangerous accessor. use it only if you know what you're doing.
	 */
	DATA_TYPE * getRawList() { return m_data; }

	/**same as above, but const correct. less dangerous. */
	const DATA_TYPE * getRawListConst() const { return m_data; }

	/**
	 * @param a_index is overwritten by the next element, which is 
	 * overwritten by the next element, and so on, till the last element
	 * @note this operation is memory intensive!
	 */
	void remove(int const a_index) {
		moveDown(a_index, -1, size());
		allocateToSize(m_allocated-1);
	}

	/** 
	 * @param a_index where to insert a_value. shifts elements in the vector.
	 * @note this operation is memory intensive!
	 */
	void insert(int const a_index, DATA_TYPE const & a_value) {
		setSize(m_size+1);
		moveUp(m_data, a_index, 1, size());
		set(a_index, a_value);
	}

	/** swaps to elements in the vector */
	void swap(int const a_index0, int const a_index1) {
		DATA_TYPE temp = get(a_index0);
		set(a_index0, get(a_index1));
		set(a_index1, temp);
	}


	/** @return index of 1st a_value range: [a_startingIndex, a_endingIndex). uses == TODO move a_list to initial param*/
	static int indexOf(DATA_TYPE const & a_value, const DATA_TYPE * const a_list, int const a_startingIndex, int const a_endingIndex) {
		for (int i = a_startingIndex; i < a_endingIndex; ++i) {
			if (a_list[i] == a_value) { return i; }
		}
		return -1;
	}

	/** @return index of 1st a_value at or after a_startingIndex. uses == */
	int indexOf(DATA_TYPE const & a_value, int const a_startingIndex, int const a_endingIndex) const {
		for(int i = a_startingIndex; i < a_endingIndex; ++i) {
			if(get(i) == a_value) { return i; }
		}
		return -1;
	}

	/** @return index of 1st a_value at or after a_startingIndex. uses operator== */
	int indexOf(DATA_TYPE const & a_value, int const a_startingIndex) const {
		return indexOf(a_value, a_startingIndex, size());
	}

	/** @return index of 1st a_value at or after a_startingIndex. uses operator== */
	int indexOf(DATA_TYPE const & a_value) const {
		return indexOf(a_value, 0, size());
	}

	/**
	 * will only work correctly if the data is sorted.
	 * @return the index of the given value, -1 if the value is not in the list
	 */
	int indexOfWithBinarySearch(DATA_TYPE const & a_value, int const a_first, int const a_limit) const {
		if(m_allocated) {
			int first = a_first, last = a_limit;
			while (first <= last) {
				int mid = (first + last) / 2;  // compute mid point.
				if (a_value > get(mid))
					first = mid + 1;  // repeat search in top half.
				else if (a_value < get(mid))
					last = mid - 1; // repeat search in bottom half.
				else return mid;     // found it. return position
			}
		}
		return -1;    // failed to find key
	}

	void setAll(DATA_TYPE const & a_value) {
		for(int i = 0; i < size(); ++i) { set(i, a_value); }
	}

	/** @return true if the given range is in order */
	bool isSorted(const int a_startIndex, const int a_endIndex) const {
		bool sorted = true;
		for (int i = a_startIndex + 1; sorted && i < a_endIndex; ++i) {
			sorted = get(i-1) < get(i);
		}
		return sorted;
	}

	bool isSorted() const { return isSorted(0, size()); }

	void reverse() {
		int half = size() / 2;
		for (int i = 0; i < half; ++i) {
			swap(i, size() - 1 - i);
		}
	}
	static void reverse(DATA_TYPE * const list, const int count) {
		int half = count / 2;
		DATA_TYPE temp;
		for (int i = 0; i < half; ++i) {
			temp = list[i];
			list[i] = list[count - 1 - i];
			list[count - 1 - i] = temp;
		}
	}

	void quicksort() { quicksort(0, size()-1); }

	/**
	* Quicksort. requires operator<. http://www.sorting-algorithms.com/static/QuicksortIsOptimal.pdf
	* @param first - The start of the sequence to be sorted (inclusive).
	* @param last - The end of the sequence to be sorted (inclusive).
	*/
	void quicksort(const int first, const int last) {
		int i = first - 1, j = last;
		DATA_TYPE v = get(last);
		if (last <= first) return;
		do {
			while (get(++i) < v);
			while (v < get(--j)) if (j == first) break;
			if (i >= j) break;
			swap(i, j);
		} while (true);
		swap(i, last);
		quicksort(first, i - 1);
		quicksort(i + 1, last);
	}

	bool containsDuplicates() const {
		for (int i = 0; i < size(); ++i) { if (indexOf(get(i), i + 1) >= 0) { return true; } }
		return false;
	}

	/** @return true if superset has every element in subset */
	static bool containsAll(const DATA_TYPE * superset, const int supersetCount, const DATA_TYPE * subset, const int subsetCount) {
		bool found = false;
		for (int i = 0; i < subsetCount; ++i) {
			found = indexOf(subset[i], superset, 0, supersetCount) >= 0;
			if (!found) break;
		}
		return found;
	}
};
