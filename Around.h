#include<cmath>
using namespace std;

#ifndef AROUND
#define AROUND

template<class T>
class Around
{
	public:
		Around();					//Default constructor
		Around(T);					//Constructor with x+/-0
		Around(T, T);					//Constructor with x+/-y
		Around(const Around&, const Around&);		//Constructor with two Around<T>s x+/-sqrt(deltax^2+y^2), to accomadate error of error integration accumalation
		Around(const Around&);				//Copy constructor
		Around<T>& operator=(Around<T>);		//Assignment
		Around<T>& operator+=(Around<T>);		//Add and assign
		Around<T>& operator-=(Around<T>);		//Subtract and assign
		Around<T> operator+(const Around<T>) const;	//Summation
		Around<T> operator+(T);
		Around<T> operator-(Around<T>);			//Difference
		Around<T> operator-(T);
		Around<T> operator*(Around<T>);			//Product
		Around<T> operator*(T);
		Around<T> operator/(Around<T>) const;		//Quotient with uncertain divisor
		Around<T> operator/(Around<T>);			//Quotient with uncertain divisor
		Around<T> operator/(T);				//Quotient with exact divisor
		Around<T> abs(Around&);				//Absolute value
		Around<T> abs();
		bool operator==(Around<T>);			//Equality test
		bool operator>=(Around<T>);			//Greater than or equal test
		bool operator>(const Around<T>) const;		//Greater than test
		bool operator<(const Around<T>) const;		//Less than test
		bool operator==(T);				//Equality to exact value. I wanted to be clever with if the value is within a standard deviation, call it equal. Caused problems in the algorithm cutting off early because small value appeared to be zero when more distance required
		bool operator>=(T);
		bool operator>(T);
		bool isnan();					//Is either part of Around<T>() nan?
		T Value();					//Return value
		T Error();					//Return error estimate aka absolute error
		T RelErr();					//Return the relative error (error/value)
		ostream& operator<<(ostream&, const Around<T>&);//Write Around<T> to stream conformal to Mathematica standard
	private:
		T value;
		T error;
};

template<class T>
ostream& operator<<(ostream& os, const Around<T>& A)
{
	//os << "Around<T>[" << A.value << "," << A.error << "]" << flush;
	os << A.value << "," << A.error << flush;
	return(os);
}

template<class T>
Around<T>::Around()
{
	value = 0;
	error = 0;
}

template<class T>
Around<T>::Around(T V)
{
	value = V;
	error = 0;
}

template<class T>
Around<T>::Around(T V, T E)
{
	value = V;
	error = E;
}

template<class T>
Around<T>::Around(const Around& V, const Around& E)
{
	value = V.value;
	error = sqrt(pow(V.error,2)+pow(E.value,2));
}

template<class T>
Around<T>::Around(const Around& A)
{
	value = A.value;
	error = A.error;
}

template<class T>
Around<T>& Around<T>::operator=(Around<T> A)
{
	error = A.error;
	value = A.value;
	return(*this);
}

template<class T>
Around<T>& Around<T>::operator+=(Around<T> A)
{
	error = sqrt(pow(error,2)+pow(A.error,2));
	value += A.value;
	return(*this);
}

template<class T>
Around<T>& Around<T>::operator-=(Around<T> A)
{
	error = sqrt(pow(error,2)+pow(A.error,2));
	value -= A.value;
	return(*this);
}

template<class T>
Around<T> Around<T>::operator+(const Around<T> A) const
{
	Around<T> B;
	B.error = sqrt(pow(error,2)+pow(A.error,2));
	B.value = value + A.value;
	return(B);
}

template<class T>
Around<T> Around<T>::operator+(T A)
{
	Around<T> B;
	B.error = error;
	B.value = value + A;
	return(B);
}

template<class T>
Around<T> operator+(T A, Around<T> B)	//exact calling sum, turn it Around<T> have the uncertain call the sum
{
	return(Around<T>(A,0)+B);
}

template<class T>
Around<T> Around<T>::operator-(Around<T> A)
{
	Around<T> B;
	B.error = sqrt(pow(error,2)+pow(A.error,2));
	B.value = value - A.value;
	return(B);
}

template<class T>
Around<T> Around<T>::operator-(T A)
{
	Around<T> B;
	B.error = error;
	B.value = value - A;
	return(B);
}

template<class T>
Around<T> operator-(T A, Around<T> B)
{
	return(Around<T>(A,0)-B);
}

template<class T>
Around<T> Around<T>::operator*(Around<T> A)
{
	Around<T> B;
	B.error = sqrt(pow(A.value*error,2)+pow(value*A.error,2));
	B.value = value * A.value;
	return(B);
}

template<class T>
Around<T> Around<T>::operator*(T A)
{
	Around<T> B;
	B.error = error * A;
	B.value = value * A;
	return(B);
}

template<class T>
Around<T> operator*(T A, Around<T> B)
{
	return(B*A);
}

template<class T>
Around<T> Around<T>::operator/(Around<T> A)
{
	Around<T> B;
	B.error = sqrt(pow(error/A.value,2)+pow(value*A.error/pow(A.value,2),2));
	B.value = value / A.value;
	return(B);
}

template<class T>
Around<T> Around<T>::operator/(Around<T> A) const
{
	Around<T> B;
	B.error = sqrt(pow(error/A.value,2)+pow(value*A.error/pow(A.value,2),2));
	B.value = value / A.value;
	return(B);
}

template<class T>
Around<T> Around<T>::operator/(T A)
{
	Around<T> B;
	B.error = error / A;
	B.value = value / A;
	return(B);
}

template<class T>
Around<T> operator/(T A, Around<T> B)	//Quotent with exact dividend.
{
	return(Around<T>(A,0)/B);
}

template<class T>
Around<T> abs(Around<T> A)
{
	return(A.abs());
}

template<class T>
Around<T> Around<T>::abs(Around& A)
{
	return(A.abs());
}

template<class T>
Around<T> Around<T>::abs()
{
	return(Around<T>(std::abs(value), error));
}

template<class T>
bool Around<T>::operator==(Around<T> A)
{
	if(std::abs(A.value-this->value) < sqrt(pow(this->error,2)+pow(A.error,2)))
		return(true);
	return(false);
}

template<class T>
bool Around<T>::operator>(const Around<T> A) const
{
	if(value > A.value && std::abs(A.value-this->value) > sqrt(pow(this->error,2)+pow(A.error,2)))
		return(true);
	return(false);
}

template<class T>
bool Around<T>::operator<(const Around<T> A) const
{
	if(value < A.value && std::abs(A.value-this->value) > sqrt(pow(this->error,2)+pow(A.error,2)))
		return(true);
	return(false);
}

template<class T>
bool Around<T>::operator>=(Around<T> A)
{
	if(*this > A || *this == A)
		return(true);
	return(false);
}

template<class T>
bool Around<T>::operator==(T A)
{
	if(value == A)
		return(true);
	return(false);
}

template<class T>
bool Around<T>::operator>(T A)
{
	if(value > A)
		return(true);
	return(false);
}

template<class T>
bool Around<T>::operator>=(T A)
{
	if(*this > A || *this == A)
		return(true);
	return(false);
}

template<class T>
bool isnan(Around<T> A)
{
	return(A.isnan());
}

template<class T>
bool Around<T>::isnan()
{
	return(std::isnan(value) || std::isnan(error));
}

template<class T>
T Around<T>::Value()
{
	return(value);
}

template<class T>
T Around<T>::Error()
{
	return(error);
}

template<class T>
T Around<T>::RelErr()
{
	return(std::abs(error/value));
}

#endif
