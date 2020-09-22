#pragma once

template <class T>

class DynamicArray
{
public:
	DynamicArray<T>(size_t initialCapacity) {
		data = new T[initialCapacity];
		capacity = initialCapacity;
		size = 0;
	}

	void setCapacity(size_t capacity) {
		delete[] data;
		data = new T[capacity];
		this->capacity = capacity;
	}

	void insertLast(T element) {		
		if (size == capacity)
			setCapacity(capacity * 2);
		data[size++] = element;
	}

	void insertLast_unsafe(T element) {
		data[size++] = element;
	}

	void clear() {
		size = 0;
	}

	size_t size() {
		return size;
	}

	size_t capacity() {
		return capacity;
	}

	~DynamicArray()
	{
		delete[] data;
	}

	T& operator[](size_t idx) {
		return data[idx];
	}

	const T operator[](size_t idx) {
		return data[idx];
	}

private:
	size_t capacity;
	size_t size;
	T *data;
};

