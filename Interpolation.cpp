#include<iostream>
#include<cmath>
#include"Interpolation.h"
using namespace std;

long double Basis0(long double x)	//Technically, Basis0-3 should only be called if x>=0, but if you're doing extrapolation, this is how it is done.
{
	if(0 <= x && x <= 2)
		return(-pow((x-2.)/2.,3));
	return(0);
}

long double Basis1(long double x)
{
	if(0 <= x && x < 2)
		return(x*(19.*pow(x,2)-90.*x+108.)/72.);
	else if(2 <= x && x <= 3)
		return(-pow(x-3,3)/9.);
	return(0);
}

long double Basis2(long double x)
{
	if(0 <= x && x < 2)
		return(-pow(x,2)*(13.*x-36.)/72.);
	else if(2 <= x && x < 3)
		return(23.*pow(x,3)/72.-2.5*pow(x,2)+6.*x-4.);
	else if(3 <= x && x <= 4)
		return(-pow((x-4.)/2.,3));
	return(0);
}

long double Basis3(long double x)
{
	if(0 <= x && x < 2)
		return(pow(x,3)/24.);
	else if(2 <= x && x < 3)
		return(-3.*pow(x/2.,3)+2.5*pow(x,2)-5.*x+10./3.);
	else if(3 <= x && x < 4)
		return(11.*pow(x,3)/24.-5*pow(x,2)+17.5*x-115./6.);
	else if(4 <= x && x <= 5)
		return(-pow((x-5.),3)/6.);
	return(0);
}

long double Basisn(long double x)
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

long double Basis_Matrix(int i, int j, int length) //i is the row from top, j is the column, and length is the length of the interpolation.
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

void Matrix_Augment_Solve_Column(long double** Data_Ctrl, long double Matrix[][4], int length)
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
	for(j = 0; j <= length; j++)
	{
		Data_Ctrl[1][j] = normalize*(multiple*Data_Ctrl[0][j]+Data_Ctrl[1][j]);
		Data_Ctrl[length-1][j] = normalize*(multiple*Data_Ctrl[length][j]+Data_Ctrl[length-1][j]);
	}
	Matrix[1][0] = Matrix[length-1][3] = 0;
	Matrix[1][1] = Matrix[length-1][2] = 1;
	Matrix[1][2] = Matrix[length-1][1] = 23./37.;
	Matrix[1][3] = Matrix[length-1][0] = 3./37.;

	//First and last 4 rows have to work the same way every time
	//Apply Row 1 and length-1 to Row 2 and length-2
	multiple = -1./9.;
	normalize = 37./18.;
	for(j = 0; j <= length; j++)
	{
		Data_Ctrl[2][j] = normalize*(multiple*Data_Ctrl[1][j]+Data_Ctrl[2][j]);
		Data_Ctrl[length-2][j] = normalize*(multiple*Data_Ctrl[length-1][j]+Data_Ctrl[length-2][j]);
	}
	Matrix[2][0] = Matrix[length-2][2] = 0;
	Matrix[2][1] = Matrix[length-2][1] = 1;
	Matrix[2][2] = Matrix[length-2][0] = 2./3.;
	Matrix[2][3] = Matrix[length-2][3] = 0;

	//First and last 4 rows have to work the same way every time
	//Apply Row 2 and length-2 to Row 3 and length-3
	multiple = -1./8.;
	normalize = 8./5.;
	for(j = 0; j <= length; j++)
	{
		Data_Ctrl[3][j] = normalize*(multiple*Data_Ctrl[2][j]+Data_Ctrl[3][j]);
		Data_Ctrl[length-3][j] = normalize*(multiple*Data_Ctrl[length-2][j]+Data_Ctrl[length-3][j]);
	}
	Matrix[3][0] = Matrix[length-3][2] = 0;
	Matrix[3][1] = Matrix[length-3][1] = 1;
	Matrix[3][2] = Matrix[length-3][0] = 4./15.;
	Matrix[3][3] = Matrix[length-3][3] = 0;

	//Interior rows going to center
	multiple = -1./6.;	//The factor will always be -1/6
	for(i = 3; i < float(length)/2.-1; i++)
	{
		normalize = 1./(Matrix[i][2]*multiple+Matrix[i+1][1]);
		for(j = 0; j <= length; j++)
		{
			Data_Ctrl[i+1][j] = normalize*(multiple*Data_Ctrl[i][j]+Data_Ctrl[i+1][j]);
			Data_Ctrl[length-i-1][j] = normalize*(multiple*Data_Ctrl[length-i][j]+Data_Ctrl[length-i-1][j]);
		}
		Matrix[i+1][2] = Matrix[length-i-1][0] = Matrix[i+1][2]*normalize;
		Matrix[i+1][0] = Matrix[length-i-1][2] = 0;
		Matrix[i+1][1] = Matrix[length-i-1][1] = 1;
		Matrix[i+1][3] = Matrix[length-i-1][3] = 0;
	}

	//Middle one or two lines
	if(length % 2 == 0)	//Even length has one line
	{
		normalize = 1./(Matrix[i][2]*multiple*2.+Matrix[i+1][1]);
		for(j = 0; j <= length; j++)
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
		for(j = 0; j <= length; j++)
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
		for(j = 0; j <= length; j++)
		{
			Data_Ctrl[i][j] = multiple*Data_Ctrl[i+1][j]+Data_Ctrl[i][j];
			Data_Ctrl[length-i][j] = multiple*Data_Ctrl[length-i-1][j]+Data_Ctrl[length-i][j];
		}
		Matrix[i][2] = Matrix[length-i][0] = 0;
	}

	//First and last 4 rows have to work the same way every time
	//Apply Row 4 and length-4 to Row 3 and length-3
	//no-op, already done by loop

	//First and last 4 rows have to work the same way every time
	//Apply Row 3 and length-3 to Row 2 and length-2
	multiple = -2./3.;
	for(j = 0; j <= length; j++)
	{
		Data_Ctrl[2][j] = multiple*Data_Ctrl[3][j]+Data_Ctrl[2][j];
		Data_Ctrl[length-2][j] = multiple*Data_Ctrl[length-3][j]+Data_Ctrl[length-2][j];
	}
	Matrix[2][2] = Matrix[length-2][0] = 0;

	//First and last 4 rows have to work the same way every time
	//Apply Row 3 and length-3 to Row 1 and length-1
	multiple = -3./37.;
	for(j = 0; j <= length; j++)
	{
		Data_Ctrl[1][j] = multiple*Data_Ctrl[3][j]+Data_Ctrl[1][j];
		Data_Ctrl[length-1][j] = multiple*Data_Ctrl[length-3][j]+Data_Ctrl[length-1][j];
	}
	Matrix[1][3] = Matrix[length-1][0] = 0;

	//First and last 4 rows have to work the same way every time
	//Apply Row 2 and length-2 to Row 1 and length-1
	multiple = -23./37.;
	for(j = 0; j <= length; j++)
	{
		Data_Ctrl[1][j] = multiple*Data_Ctrl[2][j]+Data_Ctrl[1][j];
		Data_Ctrl[length-1][j] = multiple*Data_Ctrl[length-2][j]+Data_Ctrl[length-1][j];
	}
	Matrix[1][2] = Matrix[length-1][1] = 0;

	//First and last 4 rows have to work the same way every time
	//Row 0 and length are alrady in their finished state
	//no-op
}

void Matrix_Augment_Solve_Row(long double** Data_Ctrl, long double Matrix[][4], int length)
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
	for(j = 0; j <= length; j++)
	{
		Data_Ctrl[j][1] = normalize*(multiple*Data_Ctrl[j][0]+Data_Ctrl[j][1]);
		Data_Ctrl[j][length-1] = normalize*(multiple*Data_Ctrl[j][length]+Data_Ctrl[j][length-1]);
	}
	Matrix[1][0] = Matrix[length-1][3] = 0;
	Matrix[1][1] = Matrix[length-1][2] = 1;
	Matrix[1][2] = Matrix[length-1][1] = 23./37.;
	Matrix[1][3] = Matrix[length-1][0] = 3./37.;

	//First and last 4 rows have to work the same way every time
	//Apply Row 1 and length-1 to Row 2 and length-2
	multiple = -1./9.;
	normalize = 37./18.;
	for(j = 0; j <= length; j++)
	{
		Data_Ctrl[j][2] = normalize*(multiple*Data_Ctrl[j][1]+Data_Ctrl[j][2]);
		Data_Ctrl[j][length-2] = normalize*(multiple*Data_Ctrl[j][length-1]+Data_Ctrl[j][length-2]);
	}
	Matrix[2][0] = Matrix[length-2][2] = 0;
	Matrix[2][1] = Matrix[length-2][1] = 1;
	Matrix[2][2] = Matrix[length-2][0] = 2./3.;
	Matrix[2][3] = Matrix[length-2][3] = 0;

	//First and last 4 rows have to work the same way every time
	//Apply Row 2 and length-2 to Row 3 and length-3
	multiple = -1./8.;
	normalize = 8./5.;
	for(j = 0; j <= length; j++)
	{
		Data_Ctrl[j][3] = normalize*(multiple*Data_Ctrl[j][2]+Data_Ctrl[j][3]);
		Data_Ctrl[j][length-3] = normalize*(multiple*Data_Ctrl[j][length-2]+Data_Ctrl[j][length-3]);
	}
	Matrix[3][0] = Matrix[length-3][2] = 0;
	Matrix[3][1] = Matrix[length-3][1] = 1;
	Matrix[3][2] = Matrix[length-3][0] = 4./15.;
	Matrix[3][3] = Matrix[length-3][3] = 0;

	//Interior rows going to center
	multiple = -1./6.;	//The factor will always be -1/6
	for(i = 3; i < float(length)/2.-1; i++)
	{
		normalize = 1./(Matrix[i][2]*multiple+Matrix[i+1][1]);
		for(j = 0; j <= length; j++)
		{
			Data_Ctrl[j][i+1] = normalize*(multiple*Data_Ctrl[j][i]+Data_Ctrl[j][i+1]);
			Data_Ctrl[j][length-i-1] = normalize*(multiple*Data_Ctrl[j][length-i]+Data_Ctrl[j][length-i-1]);
		}
		Matrix[i+1][2] = Matrix[length-i-1][0] = Matrix[i+1][2]*normalize;
		Matrix[i+1][0] = Matrix[length-i-1][2] = 0;
		Matrix[i+1][1] = Matrix[length-i-1][1] = 1;
		Matrix[i+1][3] = Matrix[length-i-1][3] = 0;
	}

	//Middle one or two lines
	if(length % 2 == 0)	//Even length has one line
	{
		normalize = 1./(Matrix[i][2]*multiple*2.+Matrix[i+1][1]);
		for(j = 0; j <= length; j++)
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
		for(j = 0; j <= length; j++)
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
		for(j = 0; j <= length; j++)
		{
			Data_Ctrl[j][i] = multiple*Data_Ctrl[j][i+1]+Data_Ctrl[j][i];
			Data_Ctrl[j][length-i] = multiple*Data_Ctrl[j][length-i-1]+Data_Ctrl[j][length-i];
		}
		Matrix[i][2] = Matrix[length-i][0] = 0;
	}

	//First and last 4 rows have to work the same way every time
	//Apply Row 4 and length-4 to Row 3 and length-3
	//no-op, already done by loop

	//First and last 4 rows have to work the same way every time
	//Apply Row 3 and length-3 to Row 2 and length-2
	multiple = -2./3.;
	for(j = 0; j <= length; j++)
	{
		Data_Ctrl[j][2] = multiple*Data_Ctrl[j][3]+Data_Ctrl[j][2];
		Data_Ctrl[j][length-2] = multiple*Data_Ctrl[j][length-3]+Data_Ctrl[j][length-2];
	}
	Matrix[2][2] = Matrix[length-2][0] = 0;

	//First and last 4 rows have to work the same way every time
	//Apply Row 3 and length-3 to Row 1 and length-1
	multiple = -3./37.;
	for(j = 0; j <= length; j++)
	{
		Data_Ctrl[j][1] = multiple*Data_Ctrl[j][3]+Data_Ctrl[j][1];
		Data_Ctrl[j][length-1] = multiple*Data_Ctrl[j][length-3]+Data_Ctrl[j][length-1];
	}
	Matrix[1][3] = Matrix[length-1][0] = 0;

	//First and last 4 rows have to work the same way every time
	//Apply Row 2 and length-2 to Row 1 and length-1
	multiple = -23./37.;
	for(j = 0; j <= length; j++)
	{
		Data_Ctrl[j][1] = multiple*Data_Ctrl[j][2]+Data_Ctrl[j][1];
		Data_Ctrl[j][length-1] = multiple*Data_Ctrl[j][length-2]+Data_Ctrl[j][length-1];
	}
	Matrix[1][2] = Matrix[length-1][1] = 0;

	//First and last 4 rows have to work the same way every time
	//Row 0 and length are alrady in their finished state
	//no-op
}

int main()
{
	const int Length = 10;
	//An 11x11 set of random data. The random data is uniform between -10 and 10.
	long double Data[Length+1][Length+1] = {{4.897923049837047, 0.5703313555402829, 6.230666205051591, 0.2539576271597852, -7.109081637521637, 1.9776254003749294, 8.527654870116066, -7.856180435031348, -4.633206932019711, 6.336658363126816, 5.835085557195562}, 
{4.148342602801726, -7.32922165263858, -8.67361733701042, -9.987510739149783, 9.066516015107602, -3.9928446456146496, -5.41258817158392, -8.661407424766097, 5.769180897448464, 0.225856030286522, 1.6485697780511828}, 
{-6.999119197698942, 7.579851496174911, 5.107771380252657, 4.974016931384249, 9.073370482796761, 1.7618738292628393, -1.5918164841612423, 2.5973755207389004, 4.849875899351488, 9.457495347162883, -0.2107819882079376}, 
{5.3760906136251805, -9.512279645391125, 2.011760453690087, 1.2086496057571807, 0.0999354920707951, 2.710735148381332, -4.001887178219889, 9.582253148924138, 9.218956943852554, 5.674792113762464, -1.7162612350260034}, 
{2.016068324212199, 7.9865972079525065, 5.234582512249588, 0.45821071644714095, -2.7805497064411533, -4.793108344931348, 4.433571901863672, 5.539773676114741, 9.488643398684296, 0.7059980446179779, -5.0688216502151064}, 
{-0.7204883186619213, 9.899077227546549, 6.547007863305815, -1.3729503106282728, 8.562867917810557, -2.041525348645571, 0.6815664483102033, -2.7604504695623575, 3.992996403574722, 6.552170070466502, 6.823385451146983}, 
{-2.864190755922369, 1.5805132151015044, 3.072727870206357, 4.87831384633288, -3.158323675890582, 1.0299688847879267, -9.049531513232296, 2.9086314606227077, 2.4359354387514784, 7.7227895898817565, 0.3809868849361351}, 
{-9.58819091317666, -3.604528014919417, 9.513387995162589, -5.659153659150363, 6.430823951310465, -5.3155797928006585, -4.024276467249571, -5.283238248780755, -7.451152013148349, 0.2765082145073876, 7.649097938140045}, 
{-1.9391053020599287, 9.186821373706174, -3.731911448692724, 4.7429131535093845, 0.5012370927405065, -6.776715011014559, 5.408374584793236, 1.1149118977580912, -8.97830608977553, -7.313552370534666, -4.8791759445697345}, 
{-9.335192087375844, 6.95704942027335, -5.944967176325946, -5.616101620098618, 6.593893841329198, -1.417125355092228, 2.0202388253607424, 9.590001514904344, -4.424511213736217, 2.5573496537675737, -1.8385813544547247}, 
{-0.3003219846464731, -1.2842331883503988, -3.505701528003872, 3.539982233724068, -1.5114952287728975, 9.8977341092844, -0.2612296879202347, -3.6735095071255515, 3.1745518762762615, -6.082055994318743, -2.8142696971582026}};
	//An 11x11 set of control points derived from the random data. These are the answers I'm expecting to be produced for the above data.
	long double Control[Length+1][Length+1] = {{4.897923049837047, -7.7667949466890756, 12.250620688094962, 0.8632291172261964, -11.332943501323859, 1.8140550629394192, 15.94247565181576, -14.418028449506048, -2.4039164688042054, 13.574805304361773, 5.835085557195562}, 
{20.69834394897445, -60.82259168596026, -9.934862241318342, -47.082537433232105, 32.45668039041305, -20.341750929180193, -10.374528391901812, -29.378434323330076, 48.13408051258559, -50.3511613471563, 2.4698806102807027}, 
{-23.752719412811047, 68.48742090944994, 3.977197089708101, 17.28081025342062, 18.379977954078832, -3.902070641458638, 8.751238864509844, 0.6240815477756986, -6.981761507779547, 35.38418835915953, -1.1576338928504875}, 
{11.6910601119301, -53.86954119124129, 10.474272419483942, -3.4567390352314615, -3.2777074366878702, 14.256340692115783, -21.47719958689781, 20.527583273351553, 5.004081097052377, 6.369472122598355, 0.4737503200334316}, 
{0.38407776565644647, 27.776881318352338, 1.6493858332507378, 5.963353341304135, -11.19378268683364, -12.941814068693576, 14.755866698041112, 3.335900211372523, 21.267119204869402, -15.414007356756604, -11.442780850660247}, 
{-1.1309612292826978, 5.3268170981930965, 27.72496236826563, -21.690419704200437, 36.449258225970375, -14.88080864279941, 11.072309237088541, -16.345177195712168, 11.697770716103683, 6.243636900823179, 14.884443181316925}, 
{-0.18316276049718422, 16.775571993099227, -24.17038388757557, 29.39314074903151, -31.841564876237957, 21.086737044297916, -29.788455156072615, 20.932918344946142, -6.601311311704242, 22.112482438189446, -7.154679167725529}, 
{-15.321532264262776, -51.858480761843964, 65.12203038412873, -41.03901173461332, 36.32889772681876, -19.656509215114742, 0.5099735146199961, -12.693109489384998, -12.178689465849374, 7.197606937750048, 16.020194799202}, 
{10.360705872738679, 84.22986617181519, -77.3084874896572, 51.62537892677248, -27.27738129195919, -10.834858239969604, 21.97368517867889, 4.375507071448398, -15.871570789781398, -19.398563226767966, -20.04874813339025}, 
{-23.290880289444427, -5.556196099335599, 27.38610617929206, -49.46029917511592, 42.77457257294959, -6.8363755761665335, -12.794538640215611, 33.818983922951205, -35.13491267923053, 50.449804741085366, 8.270572768217662}, 
{-0.3003219846464731, 4.849452602101386, -12.932826201279196, 9.421121550753252, -9.100253537397531, 17.91092122619949, -3.1570267116940123, -6.850192506944852, 13.639050431275091, -19.07370774905455, -2.8142696971582026}};
	long double** Control_ptr = new long double*[Length+1];
	for(int i = 0; i < Length+1; i++)
	{
		Control_ptr[i] = new long double[Length+1];
		for(int j = 0; j < Length+1; j++)
			Control_ptr[i][j] = Control[i][j];
	}
	long double** Data_ptr = new long double*[Length+1];
	for(int i = 0; i < Length+1; i++)
	{
		Data_ptr[i] = new long double[Length+1];
		for(int j = 0; j < Length+1; j++)
			Data_ptr[i][j] = Data[i][j];
	}

	long double Matrix[Length+1][4];
	for(int i = 0; i <= Length; i++)
		for(int j = 0; j < 4; j++)
			Matrix[i][j] = Basis_Matrix(i, j, Length);
	Matrix_Augment_Solve_Column(Data_ptr, Matrix, Length);

	for(int i = 0; i <= Length; i++)
		for(int j = 0; j < 4; j++)
			Matrix[i][j] = Basis_Matrix(i, j, Length);
	Matrix_Augment_Solve_Row(Data_ptr, Matrix, Length);

	cout << endl;
	for(int i = 0; i <= Length; i++)
	{
		for(int j = 0; j <= Length; j++)
			cout << Data_ptr[i][j]-Control_ptr[i][j] << ",";
		cout << endl;
	}

	cout << endl;
	for(int i = 0; i <= Length; i++)
	{
		for(int j = 0; j <= Length; j++)
			cout << Data_ptr[i][j]/Control_ptr[i][j]-1. << ",";
		cout << endl;
	}

	Interpolation<long double> f(Control_ptr, Length+1, Length+1);
	Interpolation<long double> g(Data_ptr, Length+1, Length+1);

	cout << endl;
	for(int i = 0; i <= Length; i++)
	{
		for(int j = 0; j <= Length; j++)
			cout << f(i,j)-g(i,j) << ",";
		cout << endl;
	}

	cout << endl;
	for(int i = 0; i <= Length; i++)
	{
		for(int j = 0; j <= Length; j++)
			cout << g(i,j)/f(i,j)-1. << ",";
		cout << endl;
	}

	return(0);
}
