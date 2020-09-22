#pragma once

template <class T>

/** A class implementing an array that can grow dynamically. Very similar to std::vector.
The main difference is the insertLast_unsafe method, which is faster than insertLast (or 
std:vector's push_back) because it doesn't check wheter the array has sufficient capacity 
before inserting a new element.
Use this to save a few CPU-cycles when you can guarantee that the capacity is sufficient.
*/
class DynamicArray
{
public:
	DynamicArray<T>(size_t initialCapacity) {
		data = new T[initialCapacity];
		capacity_ = initialCapacity;
		size_ = 0;
	}

	void changeCapacity(size_t capacity) {
		if (capacity == capacity_)
			return;
		T* newData = new T[capacity];
		for (int i = 0; i < size_; i++)
			newData[i] = data[i];
		delete[] data;
		data = newData;
		capacity_ = capacity;
	}

	void insertLast(T element) {		
		if (size_ == capacity_)
			changeCapacity(capacity_ * 2);
		data[size_++] = element;
	}

	void insertLast_unsafe(T element) {
		data[size_++] = element;
	}

	void clear() {
		size_ = 0;
	}

	size_t size() {
		return size_;
	}

	size_t capacity() {
		return capacity_;
	}

	~DynamicArray()
	{
		delete[] data;
	}


	T& operator[](size_t idx) {
		return data[idx];
	}

	const T operator[](size_t idx) const {
		return data[idx];
	}

private:
	size_t capacity_;
	size_t size_;
	T *data;
};

