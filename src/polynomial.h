#ifndef __polynomial_h_
#define __polynomial_h_

#include <functional>
#include <sstream>

// Don't need classes like Ring/Field because C++ has operator overloading

// F should be a Field
template <typename F> class Monomial {
 public:
	// Instead of storing coeff as a double, store as expression tree? TODO
	F coeff;
	std::vector<int> exponent;

	Monomial() { coeff = 0; }
	Monomial(F _coeff) {
		coeff = _coeff;
	};
	Monomial(F _coeff, int ex) {
		coeff = _coeff;
		exponent.assign(ex+1, 0);
		exponent[ex] = 1;
	}
	Monomial(F _coeff, const std::vector<int> &_e) { coeff=_coeff; exponent=_e; }
	Monomial(const Monomial<F> &from) { coeff = from.coeff; exponent = from.exponent; }
	int getExponent(int xx) const {
		if (xx >= exponent.size()) return 0;
		return exponent[xx];
	}
	Monomial<F> operator*(const Monomial<F> &mon2) const {
		Monomial<F> erg;
		erg.coeff = coeff * mon2.coeff;
		auto numExp = std::max(exponent.size(), mon2.exponent.size());
		erg.exponent.reserve(numExp);
		for (size_t ex=0; ex < numExp; ++ex) {
			erg.exponent.push_back(getExponent(ex) + mon2.getExponent(ex));
		}
		return erg;
	}
	bool isDivisibleBy(Monomial<F> &monomial) {
		for (size_t ex=0; ex < exponent.size(); ++ex)
			if (exponent[ex] > monomial.getExponent(ex)) return false;
		return true;
	}
	Monomial<F> operator-() {
		Monomial<F> erg(this);
		erg.coeff = -erg.coeff;
		return erg;
	}
	bool hasSameExponent(const Monomial<F> m2) const { // rename to operator== ?
		auto len = std::min(exponent.size(), m2.exponent.size());
		for (size_t ex=0; ex < len; ++ex)
			if (exponent[ex] != m2.exponent[ex]) return false;
		for (size_t ex=len; ex < exponent.size(); ++ex)
			if (exponent[ex]) return false;
		for (size_t ex=len; ex < m2.exponent.size(); ++ex)
			if (m2.exponent[ex]) return false;
		return true;
	}

	operator std::string() const {
		std::stringstream erg;
		erg << coeff;
		erg << "^{";
		for (size_t ex=0; ex < exponent.size(); ++ex) {
			erg << exponent[ex];
			if (ex < exponent.size()-1) erg << " ";
		}
		erg << "}";
		return erg.str();
	}
};

// orders by exponent; ignores coeff
template <typename F> bool operator<(const Monomial<F> &m1, const Monomial<F> &m2)
{
	int s1 = std::accumulate(m1.exponent.begin(), m1.exponent.end(), 0);
	int s2 = std::accumulate(m2.exponent.begin(), m2.exponent.end(), 0);
	if (s1 != s2) return s1 < s2;
	if (m1.exponent.size() > m2.exponent.size()) {
		for (size_t ex=m2.exponent.size(); ex < m1.exponent.size(); ++ex)
			if (m1.exponent[ex]) return false;
	} else if (m1.exponent.size() < m2.exponent.size()) {
		for (size_t ex=m1.exponent.size(); ex < m2.exponent.size(); ++ex)
			if (m2.exponent[ex]) return true;
	}
	for (size_t ex = std::min(m1.exponent.size(), m2.exponent.size()) - 1; ex >= 0; --ex) {
		if (m1.exponent[ex] == m2.exponent[ex]) continue;
		return m1.exponent[ex] < m2.exponent[ex];
	}
	return false;
};

template <typename F> bool operator>(const Monomial<F> &m1, const Monomial<F> &m2)
{ return m2 < m1; }

// F should be a Field
template <typename F> class Polynomial {
 public:
	std::set< Monomial<F> > monomials;
	Polynomial() {}
	Polynomial(F coeff) {
		if (coeff != 0) addMonomial(Monomial<F>(coeff));
	}
	Polynomial(F coeff, int ex) {
		if (coeff != 0) addMonomial(Monomial<F>(coeff, ex));
	}
	Polynomial(const Monomial<F> &m1) {
		if (m1.coeff == 0) return;
		addMonomial(m1);
	}
	Polynomial(const Polynomial<F> &other) {
		for (auto &m1 : other.monomials) addMonomial(m1);
	}
	//Polynomial<F> zero() { return Polynomial<F>(F::zero()); }
	//Polynomial<F> one() { return Polynomial<F>(F::one()); }
	bool isZero() const { return monomials.size() == 0; }
	void addMonomial(Monomial<F> m)
	{
		if (m.coeff == 0) return;
		auto iter = monomials.lower_bound(m);
		if (iter != monomials.end() && iter->hasSameExponent(m)) {
			// can't modify in-place because std::set doesn't know that order won't change
			m.coeff = iter->coeff + m.coeff;
			monomials.erase(iter);
		}
		monomials.insert(m);
	}
	void addMonomial(F _coeff, int ex) {
		addMonomial(Monomial<F>(_coeff, ex));
	}
	Polynomial<F> monomialMultiply(const Monomial<F> &monom) const {
		Polynomial<F> erg;
		for (auto &m1 : monomials) {
			erg.addMonomial(monom * m1);
		}
		return erg;
	}
	void operator*=(const Polynomial<F> &arg) {
		if (isZero() || arg.isZero()) { monomials.clear(); return; }
		auto orig(*this);
		monomials.clear();
		for (auto &m1 : arg.monomials) *this += orig.monomialMultiply(m1);
	}
	Polynomial<F> operator*(const Polynomial<F> &arg) const {
		Polynomial<F> erg(*this);
		erg *= arg;
		return erg;
	}
	void operator+=(const Polynomial<F> &other) {
		auto it1 = monomials.rbegin();
		auto it2 = other.monomials.rbegin();
		Polynomial<F> erg;
		while (it1 != monomials.rend() && it2 != other.monomials.rend()) {
			auto &top1 = *it1;
			auto &top2 = *it2;
			if (top1.hasSameExponent(top2)) {
				F newCoeff = top1.coeff + top2.coeff;
				if (newCoeff != 0) erg.addMonomial(Monomial<F>(newCoeff, top1.exponent));
				++it1; ++it2;
			} else if (top1 > top2) {
				auto result = erg.monomials.insert(top1);
				if (!result.second) mxThrow("already exists in set?");
				++it1;
			} else {
				auto result = erg.monomials.insert(top2);
				if (!result.second) mxThrow("already exists in set?");
				++it2;
			}
		}
		while (it1 != monomials.rend()) {
			auto result = erg.monomials.insert(*it1);
			if (!result.second) mxThrow("already exists in set?");
			++it1;
		}
		while (it2 != other.monomials.rend()) {
			auto result = erg.monomials.insert(*it2);
			if (!result.second) mxThrow("already exists in set?");
			++it2;
		}
		monomials = erg.monomials;
	}
	operator std::string() const {
		std::stringstream erg;
		for (auto &m1 : monomials) {
			erg << std::string(m1);
			erg << " ";
		}
		return erg.str();
	};
};

#endif // __polynomial_h_
