/** \file DynamicArray.h
*	\brief Dynamic Array class similar to std::vector.
*	\author Oyvind L Rortveit
*	\date 2020
*/

#pragma once

template <class T>

/** A class implementing an array that can grow dynamically. Very similar to std::vector.
The main difference is the insertLast_unsafe method, which is faster than insertLast (or 
std:vector's push_back) because it doesn't check wheter the array has sufficient capacity 
before inserting a new element.
Use this to save a few CPU-cycles when you can guarantee that the capacity is sufficient.
Note: The clear() function does not call destructors, if destruction is necessary,
this must be taken care of by the caller.
*/
class DynamicArray
{
public:
	DynamicArray<T>(std::size_t initialCapacity) {
		data = new T[initialCapacity];
		capacity_ = initialCapacity;
		size_ = 0;
	}

	/** Change the capacity of the array. The new capacity must be
	*	at least the current size of the array, otherwise the capacity'
	*	will not be changed. Complexity is O(size()).
	*
	*	\param capacity The new capacity after call
	*/
	void changeCapacity(std::size_t capacity) {
		if (capacity <= size_)
			return;
		if (capacity == capacity_)
			return;
		T* newData = new T[capacity];
		for (int i = 0; i < size_; i++)
			newData[i] = data[i];
		delete[] data;
		data = newData;
		capacity_ = capacity;
	}

	/** Changes the capacity only if the current capacity is smaller than
	*	the new capacity.
	*/
	void setMinimumCapacity(std::size_t capacity) {
		if (capacity > capacity_)
			changeCapacity(capacity);
	}

	/** Increases the array size by one, and inserts the element
	*	at the end of  the array. Increases the capacity (doubles it)
	*	if necessary.
	*/
	void insertLast(T element) {		
		if (size_ == capacity_)
			changeCapacity(capacity_ * 2);
		data[size_++] = element;
	}

	/** Increases the array size by one, and inserts the element
	*	at the end of  the array. Does not increase the capacity
	*	if necessary, so the caller must be certain that the 
	*	current capacity is sufficient to hold the new element.
	*/
	void insertLast_unsafe(T element) {
		data[size_++] = element;
	}

	/** Sets the size of the array to zero. Does not call any destructors.
	*/
	void clear() {
		size_ = 0;
	}

	std::size_t size() const {
		return size_;
	}

	std::size_t capacity() const {
		return capacity_;
	}

	~DynamicArray()
	{
		delete[] data;
	}


	T& operator[](std::size_t idx) {
		return data[idx];
	}

	const T operator[](std::size_t idx) const {
		return data[idx];
	}

private:
	std::size_t capacity_;
	std::size_t size_;
	T *data;
};

