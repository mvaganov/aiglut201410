#pragma once

#include "templatevector.h"

/**
 * a Vector that grows in a way that is memory stable.
 * this data structure is ideal when a vector of elements is needed,
 * and the elements need to stay stationary in memory, because they 
 * are being referenced by pointers elsewhere.
 * TODO write code to force m_allocationSize to be a power of 2, so that use left-shift can replace division, and bitwise-and can replace modulo
 */
template <class DATA_TYPE>
class TemplateVectorList
{
private:
	/** a list of arrays */
	TemplateVector<DATA_TYPE*> m_allocations;
	int m_allocationSize, m_allocated, m_size;
public:
	TemplateVectorList(const int & a_allocationPageSize)
		:m_allocationSize(a_allocationPageSize), m_allocated(0), m_size(0){}
	TemplateVectorList():m_allocationSize(16),m_allocated(0),m_size(0){}
	int size() const { return m_size; }
	bool ensureCapacity(const int a_size) {
		while(a_size >= m_allocated) {
#ifdef TEMPLATEARRAY_USES_MALLOC
			DATA_TYPE* arr = (DATA_TYPE*)malloc(sizeof(DATA_TYPE)*m_allocationSize);
#else
			DATA_TYPE* arr = new DATA_TYPE[m_allocationSize];
#endif
			if(!arr) return false;
			m_allocations.add(arr);
			m_allocated += m_allocationSize;
		}
		return true;
	}
	void setSize(int const & a_size) {
		ensureCapacity(a_size);
		m_size = a_size;
	}
	void clear() { setSize(0); }
	DATA_TYPE & get(int const & a_index) {
		int arrIndex = a_index / m_allocationSize;
		int subIndex = a_index % m_allocationSize;
		return m_allocations.get(arrIndex)[subIndex];
	}
	const DATA_TYPE & getCONST(int const a_index) const {
		int arrIndex = a_index / m_allocationSize;
		int subIndex = a_index % m_allocationSize;
		return m_allocations[arrIndex][subIndex];
	}
	DATA_TYPE & operator[](int const a_index) { return get(a_index); }

	const DATA_TYPE & operator[](int const a_index) const { return getCONST(a_index); }

	DATA_TYPE & getLast() {
		return get(m_size-1);
	}
	/** cleans up memory */
	void release() {
		for(int i = 0; i < m_allocations.size(); ++i) {
#ifdef TEMPLATEARRAY_USES_MALLOC
			free(m_allocations.get(i));
#else
			delete [] m_allocations.get(i);
#endif
		}
		m_allocations.setSize(0);
		m_allocated = 0;
		m_size = 0;
	}
	~TemplateVectorList(){release();}
	void set(int const & a_index, DATA_TYPE const & a_value) {
		// complex types must overload DATA_TYPE & operator=(const DATA_TYPE &)
		get(a_index) = a_value;
	}
	void add(DATA_TYPE const & a_value) {
		ensureCapacity(m_size);
		set(m_size++, a_value);
	}
	/** adds a new empty element to the end of the vector list */
	void add() {
		ensureCapacity(m_size++);
	}
	/** adds a standard C array to this vector list */
	void add(DATA_TYPE const * const a_array, int const & a_size) {
		for(int i = 0; i < a_size; ++i) {
			add(a_array[i]);
		}
	}
	int indexOf(DATA_TYPE const & a_value, int const & a_start) const {
		int arrIndex = a_start / m_allocationSize;
		int subIndex = a_start % m_allocationSize;
		int lastList = m_size / m_allocationSize;
		int maxInList;
		bool found = false;
		for(; !found && arrIndex < m_allocations.size(); ++arrIndex)
		{
			maxInList = (arrIndex < lastList)?m_allocationSize:(a_start % m_allocationSize);
			for(; !found && subIndex < maxInList; ++subIndex)
			{
				found = m_allocations[arrIndex][subIndex] == a_value;
			}
			subIndex = 0;
		}
		return found?(arrIndex*m_allocationSize+subIndex):-1;
	}
	int indexOf(DATA_TYPE const & a_value) const {
		return indexOf(a_value, 0);
	}
	int indexOf(DATA_TYPE * const a_memoryLocation) const {
		const DATA_TYPE * start, * end;
		for(int i = 0; i < m_allocations.size(); ++i){
			start = m_allocations[i];
			end = start+m_allocationSize;
			if(a_memoryLocation >= start && a_memoryLocation < end)
			{
				int index = (size_t)a_memoryLocation - (size_t)start;
				index /= sizeof(DATA_TYPE);
				return index + i*m_allocationSize;
			}
		}
		return -1;
	}
private:
	void moveUp(int const & a_from, int const & a_offset, int const & a_last) {
		for(int i = a_last-a_offset-1; i >= a_from; --i) {
			set(i+a_offset, get(i));
		}
	}
public:
	void insert(int const & a_index, DATA_TYPE const & a_value) {
		ensureCapacity(m_size+1);
		m_size++;
		moveUp(a_index, 1, m_size);
		set(a_index, a_value);
	}
	TemplateVectorList(const TemplateVectorList<DATA_TYPE> & toCopy)
		:m_allocationSize(toCopy.m_allocationSize),m_allocated(0),m_size(0) {
		for(int i = 0; i < toCopy.size(); ++i)
			add(toCopy.getCONST(i));
	}
#ifdef CPP11_HAS_MOVE_SEMANTICS
	/** will move the data from a_array to *this */
	void moveSemantic(TemplateVectorList & a_vectorlist) {
		// overloaded move operator= should trigger here
		m_allocations.moveSemantic(a_vectorlist.m_allocations);
		m_allocationSize = a_vectorlist.m_allocationSize;
		m_allocated = a_vectorlist.m_allocated;
		m_size = a_vectorlist.m_size;
		a_vectorlist.m_allocationSize = 0;
		a_vectorlist.m_allocated = 0;
		a_vectorlist.m_size = 0;
	}

	/**
	 * move constructor, for C++11, to make the following efficient
	 * <code>TemplateVectorList<int> list(TemplateVectorList<int>());</code>
	 */
	TemplateVectorList(TemplateVectorList<DATA_TYPE> && a_vectorlist) {
		moveSemantic(a_vectorlist);
	}

	/**
	 * move assignment, for C++11, to make the following efficient
	 * <code>TemplateVectorList<int> list = TemplateVectorList<int>();</code>
	 */
	TemplateVectorList & operator=(TemplateVectorList<DATA_TYPE> && a_vectorlist) {
		release();
		moveSemantic(a_vectorlist);
		return *this;
	}
#endif

#ifdef CPP11_HAS_LAMBDA_SEMANTICS
	/** @param f execute this code for each element of this container */
	void for_each_full(std::function<void (DATA_TYPE & value, const int index)> f) {
		for(int i = 0; i < size(); ++i)
			f(TemplateVectorList<DATA_TYPE>::get(i), i);
	}
	/** @param f execute this code for each element of this container */
	void for_each(std::function<void (DATA_TYPE & value)> f) {
		for(int i = 0; i < size(); ++i)
		{
			f(TemplateVectorList<DATA_TYPE>::get(i));
		}
	}
#endif

#ifdef CPP11_HAS_INITIALIZER_LIST
	TemplateVectorList( const std::initializer_list <DATA_TYPE> & ilist )
		:m_allocated(0), m_size(0)
    {
    	m_allocationSize = ilist.size();
        setSize(m_allocationSize);
        auto it = ilist.begin();
        int index = 0;
        while( it != ilist.end() ) {
        	set(index++, *it);
            it++;
        }
    }
#endif
};
