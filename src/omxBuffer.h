#ifndef _OMXBUFFER_H_
#define _OMXBUFFER_H_

// replace with Eigen::VectorXd or similar TODO

template<typename _Tp> class omxBuffer {
	typedef _Tp value_type;
	typedef _Tp *pointer;
	typedef _Tp &reference;
	typedef const _Tp &const_reference;
	pointer _M_start;
 public:
	omxBuffer(size_t size) {
		_M_start = new _Tp[size];
	}
	omxBuffer(size_t size, const value_type &value) {
		_M_start = new _Tp[size];
		std::fill_n(_M_start, size, value);
	}
	~omxBuffer() {
		delete [] _M_start;
	}
	reference operator[](size_t __n)
	{ return *(_M_start + __n); }
	const_reference operator[](size_t __n) const
	{ return *(_M_start + __n); }
	pointer data() const
	{ return _M_start; }
};

#endif // _OMXBUFFER_H_
