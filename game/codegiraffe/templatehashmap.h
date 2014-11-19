#pragma once

#include "templatevector.h"
#include <string.h>

template <class KEY, class VALUE>
struct KeyValuePair {
	KEY k;	VALUE v;
	KeyValuePair(){}
	KeyValuePair(KEY const & k):k(k){}
	KeyValuePair(KEY const & k, VALUE const & v):k(k),v(v){}
	// operators overloaded for TemplateVector
	bool operator<(KeyValuePair<KEY,VALUE> const & kvp) const { return k < kvp.k; }
	bool operator>(KeyValuePair<KEY,VALUE> const & kvp) const { return k > kvp.k; }
	bool operator==(KeyValuePair<KEY,VALUE> const& kvp) const { return k ==kvp.k; }
	bool operator!=(KeyValuePair<KEY,VALUE> const& kvp) const { return !operator==(kvp); }
	static int hashFunction(KEY const & k) {
		// the KEY could very likely be pointer. This hash will turn memory 
		// addresses that are int-aligned into possibly consecutive keys, and
		// essentially modulo 32, accomplished here by bitwise-and-ing 32-1
		// TODO ptrdiff_t
		return (((int)k) >> (sizeof(int) >> 1)) & 31;
	}
};

template <class VALUE>
struct NameValuePair {
	const char * k;	VALUE v;
	NameValuePair():k(0){}
	NameValuePair(const char * k):k(k){}
	NameValuePair(const char * k, VALUE v):k(k),v(v){}
	// operators overloaded for TemplateVector
	bool operator<(NameValuePair<VALUE> const & kvp) const { return strcmp(k, kvp.k) < 0; }
	bool operator>(NameValuePair<VALUE> const & kvp) const { return strcmp(k, kvp.k) > 0; }
	bool operator==(NameValuePair<VALUE> const& kvp) const { return strcmp(k, kvp.k) ==0; }
	bool operator!=(NameValuePair<VALUE> const& kvp) const { return !operator==(kvp); }
	static int hashFunction(const char* const & k) { return ((k[0]) & 31); }
};

/**
 * a simple HashMap data structure, using TemplateVectors of KeyValuePairs as
 * buckets. The hash size is DEFAULT_BUCKET_SIZE, or 32. The hashFunction 
 * values between 0 and 31 by bitshifting and masking
 */
template <class KEY, class VALUE, class KVP_STRUCT>
class TemplateHashMap_BASE
{
	TemplateVector<TemplateVector<KVP_STRUCT>*> buckets;

	// using #define instead of const int to reduce templated-member ambiguities
#define __DEFAULT_BUCKET_SIZE 32

	int numElements;
public:
	int elementCount() { return numElements; }

	void clear(){
		for (int i = 0; i < buckets.size(); ++i)
			buckets.get(i).clear();
		numElements = 0;
	}
	TemplateHashMap_BASE():numElements(0) {
		buckets.setSize(__DEFAULT_BUCKET_SIZE);
		for (int i = 0; i < buckets.size(); ++i){
			TemplateVector<KVP_STRUCT> * hashBucket = new TemplateVector<KVP_STRUCT>();
			buckets.set(i, hashBucket);
			buckets.get(i)->init();
		}
	}
	/** clears the hash table (does not delete hash elements! they had better be referenced elsewhere...) */
	void release() {
		for (int i = 0; i < buckets.size(); ++i){
			TemplateVector<KVP_STRUCT>* list = buckets.get(i);
			if(list){
				delete list;
				buckets.set(i, 0);
			}
		}
		buckets.clear();
	}
	/** clears and deletes hash table elements (the elements had better be pointers!) */
	void deleteAll() {
		for (int y = 0; y < buckets.size(); ++y) {
			TemplateVector<KVP_STRUCT>* list = buckets.get(y);
			for(int x = 0; x < list->size(); ++x) {
				if(list->get(x).v){
					delete list->get(x).v;
					list->get(x).v = 0;
				}
			}
			list->clear();
			delete list;
			buckets.set(y, 0);
		}
		buckets.clear();
	}
	~TemplateHashMap_BASE(){ release(); buckets.release(); }

	/** @return the value associated with the given KEY */
	VALUE * getByKey(KEY const & k) {
		int index = KVP_STRUCT::hashFunction(k);
		TemplateVector<KVP_STRUCT> * bucket = hash.get(index);
		KVP_STRUCT kvp(k);
		int bucketIndex = bucket->indexOfWithBinarySearch(kvp);
		if(bucketIndex < 0)
			return 0;
		return &bucket->get(bucketIndex).v;
	}

	/** sets up a key/value pair association */
	void set(KEY const & k, VALUE const & v) {
		int index = KVP_STRUCT::hashFunction(k);
		TemplateVector<KVP_STRUCT> * bucket = buckets.get(index);
		int indexInserted;
		indexInserted = bucket->insertSorted(KVP_STRUCT(k,v), true);
		if(indexInserted >= 0) {
			numElements++;
		}
	}
	VALUE & operator[](KEY const & k) {
		KVP_STRUCT kvp(k);
		int index = KVP_STRUCT::hashFunction(k);
		TemplateVector<KVP_STRUCT> * bucket = buckets.get(index);
		int bucketIndex = bucket->indexOfWithBinarySearch(kvp);
		if (bucketIndex < 0) {
			bucketIndex = bucket->insertSorted(kvp, true);
			numElements++;
		}
		return bucket->operator[](bucketIndex).v;
	}

	const VALUE & operator[](KEY const & k) const {
		KVP_STRUCT kvp(k);
		int index = KVP_STRUCT::hashFunction(k);
		TemplateVector<KVP_STRUCT> * bucket = buckets.get(index);
		int bucketIndex = bucket->indexOfWithBinarySearch(kvp);
		if (bucketIndex < 0) {
			// should not be adding elements to a constant map
			int i = 0; i = 1 / i; // force a loud problem
		}
		return bucket->operator[](bucketIndex).v;
	}

	/** @return the structure that does all the work for the hash map */
	TemplateVector<TemplateVector<KVP_STRUCT>*> * getRawBuckets() {
		return buckets;
	}
#undef __DEFAULT_BUCKET_SIZE
};

template <class KEY, class VALUE>
class TemplateHashMap : public TemplateHashMap_BASE<KEY, VALUE, KeyValuePair<KEY,VALUE> >{};

template <class VALUE>
class TemplateHashMapNamed : public TemplateHashMap_BASE<const char*, VALUE, NameValuePair<VALUE> >{};
