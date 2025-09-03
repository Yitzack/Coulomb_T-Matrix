#ifndef INTERPOLATION
#define INTERPOLATION

template <class T>
class Interpolation
{
	public:
		~Interpolation();				//Deconstructor
		Interpolation();				//Default Constructor
		Interpolation(T**, int xSize, int ySize);	//Constructor with long double array
		T operator()(long double x, long double y);	//Interpolation evaluation
		void operator=(Interpolation<T>);		//Assignment operator
		int MaxX();					//Return the largest x
		int MaxY();					//Return the largest Y
		bool is_ready();
	private:
		bool ready;
		T** control_points;
		int*** offset;
		long double xRange, yRange;
		long double Basis0(long double);
		long double Basis1(long double);
		long double Basis2(long double);
		long double Basis3(long double);
		long double Basisn(long double);
		long double Basis_Wrapper(long double, int);
		long double Basis_Matrix(int, int, int);
		void Matrix_Augment_Solve_X(T**, long double[][4], int, int);
		void Matrix_Augment_Solve_Y(T**, long double[][4], int, int);
};

template <class T>
Interpolation<T>::~Interpolation()
{
	if(ready)
	{
		for(int i = 0; i <= xRange; i++)
		{
			for(int j = 0; j < yRange; j++)
			{
				delete offset[i][j];
			}
			delete offset[i];
			delete control_points[i];
		}
		delete offset;
		delete control_points;
	}
}

template <class T>
Interpolation<T>::Interpolation()
{
	ready = false;
}

template <class T>
Interpolation<T>::Interpolation(T** Data, int xSize, int ySize)	//I really wanted to derive the control points myself, but it is a lot easier to scrape them from Mathematica than figure out the matrix coefficents and inversion.
{
	xRange = xSize;
	yRange = ySize;

	control_points = new T*[xSize+1];
	offset = new int**[xSize+1];
	for(int i = 0; i < xSize+1; i++)
	{
		control_points[i] = new T[ySize+1];
		offset[i] = new int*[ySize+1];
		for(int j = 0; j < ySize+1; j++)
		{
			control_points[i][j] = Data[i][j];
			offset[i][j] = new int[2];
			offset[i][j][0] = i-5;
			offset[i][j][1] = j-5;
	}	}

	long double MatrixX[xSize+1][4];
	for(int i = 0; i <= xSize; i++)
		for(int j = 0; j < 4; j++)
			MatrixX[i][j] = Basis_Matrix(i, j, xSize);
	Matrix_Augment_Solve_X(control_points, MatrixX, xSize, ySize);

	long double MatrixY[ySize+1][4];
	for(int i = 0; i <= ySize; i++)
		for(int j = 0; j < 4; j++)
			MatrixY[i][j] = Basis_Matrix(i, j, ySize);
	Matrix_Augment_Solve_Y(control_points, MatrixY, xSize, ySize);

	ready = true;
}

template <class T>
void Interpolation<T>::operator=(Interpolation<T> A)
{
	if(A.ready)	//Assign the interpolation in A if it exists
	{
		xRange = A.xRange;
		yRange = A.yRange;

		control_points = new T*[int(xRange)+1];
		offset = new int**[int(xRange)+1];
		for(int i = 0; i <= xRange; i++)
		{
			control_points[i] = new T[int(yRange)+1];
			offset[i] = new int*[int(yRange)+1];
			for(int j = 0; j <= yRange; j++)
			{
				control_points[i][j] = A.control_points[i][j];
				offset[i][j] = new int[2];
				offset[i][j][0] = i-5;
				offset[i][j][1] = j-5;
		}	}

		ready = true;
	}
}

template <class T>
T Interpolation<T>::operator()(long double x, long double y)
{
	int Basisi[4] = {4,4,4,4};		//Basis Functions in the i direction
	int Basisj[4] = {4,4,4,4};		//Basis Functions in the j direction
	long double zx[4][4], zy[4][4];	//z from the x direction and y direction
	int offset_i = -1, offset_j = -1;	//Index offsets in calling up control points, normally -1, but can be other values on the ends
	T Answer = T(0);

	if(x < 2)	//Reassign the function pointers for the x/i direction
	{
		Basisi[0] = 0;
		Basisi[1] = 1;
		Basisi[2] = 2;
		Basisi[3] = 3;
		if(x < 1)
			offset_i = 0;
	}
	else if(x < 3)
	{
		Basisi[0] = 1;
		Basisi[1] = 2;
		Basisi[2] = 3;
	}
	else if(x < 4)
	{
		Basisi[0] = 2;
		Basisi[1] = 3;
	}
	else if(x < 5)
	{
		Basisi[0] = 3;
	}
	else if(xRange-2 <= x)
	{
		Basisi[0] = 3;
		Basisi[1] = 2;
		Basisi[2] = 1;
		Basisi[3] = 0;
	}
	else if(xRange-3 <= x)
	{
		Basisi[1] = 3;
		Basisi[2] = 2;
		Basisi[3] = 1;
	}
	else if(xRange-4 <= x)
	{
		Basisi[2] = 3;
		Basisi[3] = 2;
	}
	else if(xRange-5 <= x)
	{
		Basisi[3] = 3;
	}

	if(y < 2)	//Reassign the function pointers for the x/i direction
	{
		Basisj[0] = 0;
		Basisj[1] = 1;
		Basisj[2] = 2;
		Basisj[3] = 3;
		if(y < 1)
			offset_j = 0;
	}
	else if(y < 3)
	{
		Basisj[0] = 1;
		Basisj[1] = 2;
		Basisj[2] = 3;
	}
	else if(y < 4)
	{
		Basisj[0] = 2;
		Basisj[1] = 3;
	}
	else if(y < 5)
	{
		Basisj[0] = 3;
	}
	else if(yRange-2 <= y)
	{
		Basisj[0] = 3;
		Basisj[1] = 2;
		Basisj[2] = 1;
		Basisj[3] = 0;
	}
	else if(yRange-3 <= y)
	{
		Basisj[1] = 3;
		Basisj[2] = 2;
		Basisj[3] = 1;
	}
	else if(yRange-4 <= y)
	{
		Basisj[2] = 3;
		Basisj[3] = 2;
	}
	else if(yRange-5 <= y)
	{
		Basisj[3] = 3;
	}

	if(int(x)+3+offset_i >= xRange+1)
		offset_i = xRange-4-int(x)-offset_i;
	if(int(y)+3+offset_j >= yRange+1)
		offset_j = yRange-4-int(y)-offset_j;

	for(int i_count = 3; i_count >= 0; i_count--)	//Evaluate the Basis Functions
		for(int j_count = 3; j_count >= 0; j_count--)
		{
			if(Basisj[j_count] != 4 && y < 5)
				zy[i_count][j_count] = Basis_Wrapper(y, Basisj[j_count]);
			else if(Basisj[j_count] != 4 && yRange-5 <= y)
				zy[i_count][j_count] = Basis_Wrapper(yRange-y, Basisj[j_count]);
			else
				zy[i_count][j_count] = Basis_Wrapper(y-offset[int(x)][int(y)][1]-j_count, Basisj[j_count]);	//Some how the fractional part of the arguments were trading 

			if(Basisi[i_count] != 4 && x < 5)
				zx[i_count][j_count] = Basis_Wrapper(x, Basisi[i_count]);
			else if(Basisi[i_count] != 4 && xRange-5 <= x)
				zx[i_count][j_count] = Basis_Wrapper(xRange-x, Basisi[i_count]);
			else
				zx[i_count][j_count] = Basis_Wrapper(x-offset[int(x)][int(y)][0]-i_count, Basisi[i_count]);

			Answer += zx[i_count][j_count]*zy[i_count][j_count]*control_points[int(x)+i_count+offset_i][int(y)+j_count+offset_j];
		}

	return(Answer);
}

template <class T>
long double Interpolation<T>::Basis_Wrapper(long double x, int Basis)
{
	switch(Basis)
	{
	case 0:
		return(Basis0(x));
		break;
	case 1:
		return(Basis1(x));
		break;
	case 2:
		return(Basis2(x));
		break;
	case 3:
		return(Basis3(x));
		break;
	case 4:
	default:
		return(Basisn(x));
		break;
	}
}

template <class T>
long double Interpolation<T>::Basis0(long double x)	//Technically, Basis0-3 should only be called if x>=0, but if you're doing extrapolation, this is how it is done.
{
	if(x <= 2)
		return(-pow((x-2.)/2.,3));
	return(0);
}

template <class T>
long double Interpolation<T>::Basis1(long double x)
{
	if(x < 2)
		return(x*(19.*pow(x,2)-90.*x+108.)/72.);
	else if(x <= 3)
		return(-pow(x-3,3)/9.);
	return(0);
}

template <class T>
long double Interpolation<T>::Basis2(long double x)
{
	if(x < 2)
		return(-pow(x,2)*(13.*x-36.)/72.);
	else if(x < 3)
		return(23.*pow(x,3)/72.-2.5*pow(x,2)+6.*x-4.);
	else if(x <= 4)
		return(-pow((x-4.)/2.,3));
	return(0);
}

template <class T>
long double Interpolation<T>::Basis3(long double x)
{
	if(x < 2)
		return(pow(x,3)/24.);
	else if(x < 3)
		return(-3.*pow(x/2.,3)+2.5*pow(x,2)-5.*x+10./3.);
	else if(x < 4)
		return(11.*pow(x,3)/24.-5*pow(x,2)+17.5*x-115./6.);
	else if(x <= 5)
		return(-pow((x-5.),3)/6.);
	return(0);
}

template <class T>
long double Interpolation<T>::Basisn(long double x)
{
	if(x < 2)
		return(0);
	else if(x < 3)
		return(pow(x-2.,3)/6.);
	else if(x < 4)
		return(-pow(x,3)/2.+5*pow(x,2)-16.*x+50./3.);
	else if(x < 5)
		return(pow(x,3)/2.-7*pow(x,2)+32.*x-142./3.);
	else if(x <= 6)
		return(-pow((x-6.),3)/6.);
	return(0);
}

template <class T>
int Interpolation<T>::MaxX()
{
	return(xRange);
}

template <class T>
int Interpolation<T>::MaxY()
{
	return(yRange);
}

template <class T>
bool Interpolation<T>::is_ready()
{
	return(ready);
}

template <class T>
long double Interpolation<T>::Basis_Matrix(int i, int j, int length) //i is the row from top, j is the column, and length is the length of the interpolation.
{
	static const long double endcaps[4][5] = {
		{1,	0,	0,	0,	0},
		{1./8.,	37./72.,23./72.,1./24.,	0},
		{1./9.,	5./9.,	1./3.,	0,	0},
		{1./8.,	17./24.,1./6.,	0,	0}};

	static const long double interior[4] = {1./6., 2./3., 1./6., 0};

	if(i < 4)	//Top left corner of matrix
		return(endcaps[i][j]);
	else if(i > length - 2)	//last 2 lines of bottom right corner of matrix
			return(endcaps[length-i][3-j]);
	else if(i > length - 4)	//first 2 lines of bottom right corner of matrix
			return(endcaps[length-i][2-j]);
	return(interior[j]);	//Middle of matrix
}

template <class T>
void Interpolation<T>::Matrix_Augment_Solve_X(T** Data_Ctrl, long double Matrix[][4], int lengthX, int lengthY)
{
	long double multiple;
	long double normalize;
	int i, j;

	//First and last 4 rows have to work the same way every time
	//Normalize Row 0 and length
	//no-op Its alreay normalized and ready to go

	//First and last 4 rows have to work the same way every time
	//Apply Row 0 and length to Row 1 and length-1
	multiple = -1./8.;
	normalize = 72./37.;
	for(j = 0; j <= lengthY; j++)
	{
		Data_Ctrl[1][j] = normalize*(multiple*Data_Ctrl[0][j]+Data_Ctrl[1][j]);
		Data_Ctrl[lengthX-1][j] = normalize*(multiple*Data_Ctrl[lengthX][j]+Data_Ctrl[lengthX-1][j]);
	}
	Matrix[1][0] = Matrix[lengthX-1][3] = 0;
	Matrix[1][1] = Matrix[lengthX-1][2] = 1;
	Matrix[1][2] = Matrix[lengthX-1][1] = 23./37.;
	Matrix[1][3] = Matrix[lengthX-1][0] = 3./37.;

	//First and last 4 rows have to work the same way every time
	//Apply Row 1 and length-1 to Row 2 and length-2
	multiple = -1./9.;
	normalize = 37./18.;
	for(j = 0; j <= lengthY; j++)
	{
		Data_Ctrl[2][j] = normalize*(multiple*Data_Ctrl[1][j]+Data_Ctrl[2][j]);
		Data_Ctrl[lengthX-2][j] = normalize*(multiple*Data_Ctrl[lengthX-1][j]+Data_Ctrl[lengthX-2][j]);
	}
	Matrix[2][0] = Matrix[lengthX-2][2] = 0;
	Matrix[2][1] = Matrix[lengthX-2][1] = 1;
	Matrix[2][2] = Matrix[lengthX-2][0] = 2./3.;
	Matrix[2][3] = Matrix[lengthX-2][3] = 0;

	//First and last 4 rows have to work the same way every time
	//Apply Row 2 and length-2 to Row 3 and length-3
	multiple = -1./8.;
	normalize = 8./5.;
	for(j = 0; j <= lengthY; j++)
	{
		Data_Ctrl[3][j] = normalize*(multiple*Data_Ctrl[2][j]+Data_Ctrl[3][j]);
		Data_Ctrl[lengthX-3][j] = normalize*(multiple*Data_Ctrl[lengthX-2][j]+Data_Ctrl[lengthX-3][j]);
	}
	Matrix[3][0] = Matrix[lengthX-3][2] = 0;
	Matrix[3][1] = Matrix[lengthX-3][1] = 1;
	Matrix[3][2] = Matrix[lengthX-3][0] = 4./15.;
	Matrix[3][3] = Matrix[lengthX-3][3] = 0;

	//Interior rows going to center
	multiple = -1./6.;	//The factor will always be -1/6
	for(i = 3; i < float(lengthX)/2.-1; i++)
	{
		normalize = 1./(Matrix[i][2]*multiple+Matrix[i+1][1]);
		for(j = 0; j <= lengthY; j++)
		{
			Data_Ctrl[i+1][j] = normalize*(multiple*Data_Ctrl[i][j]+Data_Ctrl[i+1][j]);
			Data_Ctrl[lengthX-i-1][j] = normalize*(multiple*Data_Ctrl[lengthX-i][j]+Data_Ctrl[lengthX-i-1][j]);
		}
		Matrix[i+1][2] = Matrix[lengthX-i-1][0] = Matrix[i+1][2]*normalize;
		Matrix[i+1][0] = Matrix[lengthX-i-1][2] = 0;
		Matrix[i+1][1] = Matrix[lengthX-i-1][1] = 1;
		Matrix[i+1][3] = Matrix[lengthX-i-1][3] = 0;
	}

	//Middle one or two lines
	if(lengthX % 2 == 0)	//Even length has one line
	{
		normalize = 1./(Matrix[i][2]*multiple*2.+Matrix[i+1][1]);
		for(j = 0; j <= lengthY; j++)
			Data_Ctrl[i+1][j] = normalize*(multiple*Data_Ctrl[i][j]+multiple*Data_Ctrl[i+2][j]+Data_Ctrl[i+1][j]);
		Matrix[i+1][0] = 0;
		Matrix[i+1][1] = 1;
		Matrix[i+1][2] = 0;
		Matrix[i+1][3] = 0;
	}
	else	//Odd length has two lines
	{
		multiple = -Matrix[i][2];
		normalize = 1./(Matrix[i][2]*multiple+Matrix[i+1][1]);
		for(j = 0; j <= lengthY; j++)
		{
			Data_Ctrl[i+1][j] = normalize*(multiple*Data_Ctrl[i][j]+Data_Ctrl[i+1][j]);
			Data_Ctrl[i][j] = multiple*Data_Ctrl[i+1][j]+Data_Ctrl[i][j];
		}
		Matrix[i+1][0] = Matrix[i][0] = 0;
		Matrix[i+1][1] = Matrix[i][1] = 1;
		Matrix[i+1][2] = Matrix[i][2] = 0;
		Matrix[i+1][3] = Matrix[i][3] = 0;
		i--;	//Loop starts in the middle with null action without
	}

	//Interior rows going to corners
	for(i; i >= 3; i--)
	{
		multiple = -Matrix[i][2];
		for(j = 0; j <= lengthY; j++)
		{
			Data_Ctrl[i][j] = multiple*Data_Ctrl[i+1][j]+Data_Ctrl[i][j];
			Data_Ctrl[lengthX-i][j] = multiple*Data_Ctrl[lengthX-i-1][j]+Data_Ctrl[lengthX-i][j];
		}
		Matrix[i][2] = Matrix[lengthX-i][0] = 0;
	}

	//First and last 4 rows have to work the same way every time
	//Apply Row 4 and length-4 to Row 3 and length-3
	//no-op, already done by loop

	//First and last 4 rows have to work the same way every time
	//Apply Row 3 and length-3 to Row 2 and length-2
	multiple = -2./3.;
	for(j = 0; j <= lengthY; j++)
	{
		Data_Ctrl[2][j] = multiple*Data_Ctrl[3][j]+Data_Ctrl[2][j];
		Data_Ctrl[lengthX-2][j] = multiple*Data_Ctrl[lengthX-3][j]+Data_Ctrl[lengthX-2][j];
	}
	Matrix[2][2] = Matrix[lengthX-2][0] = 0;

	//First and last 4 rows have to work the same way every time
	//Apply Row 3 and length-3 to Row 1 and length-1
	multiple = -3./37.;
	for(j = 0; j <= lengthY; j++)
	{
		Data_Ctrl[1][j] = multiple*Data_Ctrl[3][j]+Data_Ctrl[1][j];
		Data_Ctrl[lengthX-1][j] = multiple*Data_Ctrl[lengthX-3][j]+Data_Ctrl[lengthX-1][j];
	}
	Matrix[1][3] = Matrix[lengthX-1][0] = 0;

	//First and last 4 rows have to work the same way every time
	//Apply Row 2 and length-2 to Row 1 and length-1
	multiple = -23./37.;
	for(j = 0; j <= lengthY; j++)
	{
		Data_Ctrl[1][j] = multiple*Data_Ctrl[2][j]+Data_Ctrl[1][j];
		Data_Ctrl[lengthX-1][j] = multiple*Data_Ctrl[lengthX-2][j]+Data_Ctrl[lengthX-1][j];
	}
	Matrix[1][2] = Matrix[lengthX-1][1] = 0;

	//First and last 4 rows have to work the same way every time
	//Row 0 and length are alrady in their finished state
	//no-op
}

template <class T>
void Interpolation<T>::Matrix_Augment_Solve_Y(T** Data_Ctrl, long double Matrix[][4], int lengthX, int lengthY)
{
	long double multiple;
	long double normalize;
	int i, j;

	//First and last 4 rows have to work the same way every time
	//Normalize Row 0 and length
	//no-op Its alreay normalized and ready to go

	//First and last 4 rows have to work the same way every time
	//Apply Row 0 and length to Row 1 and length-1
	multiple = -1./8.;
	normalize = 72./37.;
	for(j = 0; j <= lengthX; j++)
	{
		Data_Ctrl[j][1] = normalize*(multiple*Data_Ctrl[j][0]+Data_Ctrl[j][1]);
		Data_Ctrl[j][lengthY-1] = normalize*(multiple*Data_Ctrl[j][lengthY]+Data_Ctrl[j][lengthY-1]);
	}
	Matrix[1][0] = Matrix[lengthY-1][3] = 0;
	Matrix[1][1] = Matrix[lengthY-1][2] = 1;
	Matrix[1][2] = Matrix[lengthY-1][1] = 23./37.;
	Matrix[1][3] = Matrix[lengthY-1][0] = 3./37.;

	//First and last 4 rows have to work the same way every time
	//Apply Row 1 and length-1 to Row 2 and length-2
	multiple = -1./9.;
	normalize = 37./18.;
	for(j = 0; j <= lengthX; j++)
	{
		Data_Ctrl[j][2] = normalize*(multiple*Data_Ctrl[j][1]+Data_Ctrl[j][2]);
		Data_Ctrl[j][lengthY-2] = normalize*(multiple*Data_Ctrl[j][lengthY-1]+Data_Ctrl[j][lengthY-2]);
	}
	Matrix[2][0] = Matrix[lengthY-2][2] = 0;
	Matrix[2][1] = Matrix[lengthY-2][1] = 1;
	Matrix[2][2] = Matrix[lengthY-2][0] = 2./3.;
	Matrix[2][3] = Matrix[lengthY-2][3] = 0;

	//First and last 4 rows have to work the same way every time
	//Apply Row 2 and length-2 to Row 3 and length-3
	multiple = -1./8.;
	normalize = 8./5.;
	for(j = 0; j <= lengthX; j++)
	{
		Data_Ctrl[j][3] = normalize*(multiple*Data_Ctrl[j][2]+Data_Ctrl[j][3]);
		Data_Ctrl[j][lengthY-3] = normalize*(multiple*Data_Ctrl[j][lengthY-2]+Data_Ctrl[j][lengthY-3]);
	}
	Matrix[3][0] = Matrix[lengthY-3][2] = 0;
	Matrix[3][1] = Matrix[lengthY-3][1] = 1;
	Matrix[3][2] = Matrix[lengthY-3][0] = 4./15.;
	Matrix[3][3] = Matrix[lengthY-3][3] = 0;

	//Interior rows going to center
	multiple = -1./6.;	//The factor will always be -1/6
	for(i = 3; i < float(lengthY)/2.-1; i++)
	{
		normalize = 1./(Matrix[i][2]*multiple+Matrix[i+1][1]);
		for(j = 0; j <= lengthX; j++)
		{
			Data_Ctrl[j][i+1] = normalize*(multiple*Data_Ctrl[j][i]+Data_Ctrl[j][i+1]);
			Data_Ctrl[j][lengthY-i-1] = normalize*(multiple*Data_Ctrl[j][lengthY-i]+Data_Ctrl[j][lengthY-i-1]);
		}
		Matrix[i+1][2] = Matrix[lengthY-i-1][0] = Matrix[i+1][2]*normalize;
		Matrix[i+1][0] = Matrix[lengthY-i-1][2] = 0;
		Matrix[i+1][1] = Matrix[lengthY-i-1][1] = 1;
		Matrix[i+1][3] = Matrix[lengthY-i-1][3] = 0;
	}

	//Middle one or two lines
	if(lengthY % 2 == 0)	//Even length has one line
	{
		normalize = 1./(Matrix[i][2]*multiple*2.+Matrix[i+1][1]);
		for(j = 0; j <= lengthX; j++)
			Data_Ctrl[j][i+1] = normalize*(multiple*Data_Ctrl[j][i]+multiple*Data_Ctrl[j][i+2]+Data_Ctrl[j][i+1]);
		Matrix[i+1][0] = 0;
		Matrix[i+1][1] = 1;
		Matrix[i+1][2] = 0;
		Matrix[i+1][3] = 0;
	}
	else	//Odd length has two lines
	{
		multiple = -Matrix[i][2];
		normalize = 1./(Matrix[i][2]*multiple+Matrix[i+1][1]);
		for(j = 0; j <= lengthX; j++)
		{
			Data_Ctrl[j][i+1] = normalize*(multiple*Data_Ctrl[j][i]+Data_Ctrl[j][i+1]);
			Data_Ctrl[j][i] = multiple*Data_Ctrl[j][i+1]+Data_Ctrl[j][i];
		}
		Matrix[i+1][0] = Matrix[i][0] = 0;
		Matrix[i+1][1] = Matrix[i][1] = 1;
		Matrix[i+1][2] = Matrix[i][2] = 0;
		Matrix[i+1][3] = Matrix[i][3] = 0;
		i--;	//Loop starts in the middle with null action without
	}

	//Interior rows going to corners
	for(i; i >= 3; i--)
	{
		multiple = -Matrix[i][2];
		for(j = 0; j <= lengthX; j++)
		{
			Data_Ctrl[j][i] = multiple*Data_Ctrl[j][i+1]+Data_Ctrl[j][i];
			Data_Ctrl[j][lengthY-i] = multiple*Data_Ctrl[j][lengthY-i-1]+Data_Ctrl[j][lengthY-i];
		}
		Matrix[i][2] = Matrix[lengthY-i][0] = 0;
	}

	//First and last 4 rows have to work the same way every time
	//Apply Row 4 and length-4 to Row 3 and length-3
	//no-op, already done by loop

	//First and last 4 rows have to work the same way every time
	//Apply Row 3 and length-3 to Row 2 and length-2
	multiple = -2./3.;
	for(j = 0; j <= lengthX; j++)
	{
		Data_Ctrl[j][2] = multiple*Data_Ctrl[j][3]+Data_Ctrl[j][2];
		Data_Ctrl[j][lengthY-2] = multiple*Data_Ctrl[j][lengthY-3]+Data_Ctrl[j][lengthY-2];
	}
	Matrix[2][2] = Matrix[lengthY-2][0] = 0;

	//First and last 4 rows have to work the same way every time
	//Apply Row 3 and length-3 to Row 1 and length-1
	multiple = -3./37.;
	for(j = 0; j <= lengthX; j++)
	{
		Data_Ctrl[j][1] = multiple*Data_Ctrl[j][3]+Data_Ctrl[j][1];
		Data_Ctrl[j][lengthY-1] = multiple*Data_Ctrl[j][lengthY-3]+Data_Ctrl[j][lengthY-1];
	}
	Matrix[1][3] = Matrix[lengthY-1][0] = 0;

	//First and last 4 rows have to work the same way every time
	//Apply Row 2 and length-2 to Row 1 and length-1
	multiple = -23./37.;
	for(j = 0; j <= lengthX; j++)
	{
		Data_Ctrl[j][1] = multiple*Data_Ctrl[j][2]+Data_Ctrl[j][1];
		Data_Ctrl[j][lengthY-1] = multiple*Data_Ctrl[j][lengthY-2]+Data_Ctrl[j][lengthY-1];
	}
	Matrix[1][2] = Matrix[lengthY-1][1] = 0;

	//First and last 4 rows have to work the same way every time
	//Row 0 and length are alrady in their finished state
	//no-op
}

#endif
