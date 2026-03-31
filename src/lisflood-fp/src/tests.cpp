#include "lisflood.h"
#ifdef TESTING
#include "lisflood2\DataTypes.h"
void SGC2_UpdateLoadBalance(const int grid_rows, const int grid_cols_padded,
	const SubGridRowList * sub_grid_layout,
	WetDryRowBound* wet_dry_bounds);


void TestInterpolate()
{
	NUMERIC_TYPE interpData[12];
	interpData[0] = 0; //val
	interpData[1] = 0; //time

	interpData[2] = 10; //val
	interpData[3] = 10; //time

	interpData[4] = 200; //val
	interpData[5] = 20;

	interpData[6] = 10;
	interpData[7] = 30;

	interpData[8] = 10;
	interpData[9] = 40;
	interpData[10] = -1;
	interpData[11] = -1;

	NUMERIC_TYPE times[5];
	NUMERIC_TYPE values[5];
	values[0] = 0;
	times[0] = 0;

	values[1] = 10;
	times[1] = 10;

	values[2] = 200;
	times[2] = 20;

	values[3] = 10;
	times[3] = 30;

	values[4] = 10;
	times[4] = 40;


	TimeSeries ts = TimeSeries();
	ts.count = 5;
	ts.time = times;
	ts.value = values;

	NUMERIC_TYPE x, y, t;

	t = C(0.0);
	x = InterpolateTimeSeries(interpData, t);
	y = InterpolateTimeSeries(&ts, t);
	printf("X %" NUM_FMT" Y %" NUM_FMT" diff %" NUM_FMT"\n", x, y, x - y);

	t = C(3.0);
	x = InterpolateTimeSeries(interpData, t);
	y = InterpolateTimeSeries(&ts, t);
	printf("X %" NUM_FMT" Y %" NUM_FMT" diff %" NUM_FMT"\n", x, y, x - y);


	t = C(5.0);
	x = InterpolateTimeSeries(interpData, t);
	y = InterpolateTimeSeries(&ts, t);
	printf("X %" NUM_FMT" Y %" NUM_FMT" diff %" NUM_FMT"\n", x, y, x - y);

	t = C(10.0);
	x = InterpolateTimeSeries(interpData, t);
	y = InterpolateTimeSeries(&ts, t);
	printf("X %" NUM_FMT" Y %" NUM_FMT" diff %" NUM_FMT"\n", x, y, x - y);

	t = C(15.0);
	x = InterpolateTimeSeries(interpData, t);
	y = InterpolateTimeSeries(&ts, t);
	printf("X %" NUM_FMT" Y %" NUM_FMT" diff %" NUM_FMT"\n", x, y, x - y);

	t = C(20.0);
	x = InterpolateTimeSeries(interpData, t);
	y = InterpolateTimeSeries(&ts, t);
	printf("X %" NUM_FMT" Y %" NUM_FMT" diff %" NUM_FMT"\n", x, y, x - y);

	t = C(25.0);
	x = InterpolateTimeSeries(interpData, t);
	y = InterpolateTimeSeries(&ts, t);
	printf("X %" NUM_FMT" Y %" NUM_FMT" diff %" NUM_FMT"\n", x, y, x - y);

	t = C(30.0);
	x = InterpolateTimeSeries(interpData, t);
	y = InterpolateTimeSeries(&ts, t);
	printf("X %" NUM_FMT" Y %" NUM_FMT" diff %" NUM_FMT"\n", x, y, x - y);

	t = C(35.0);
	x = InterpolateTimeSeries(interpData, t);
	y = InterpolateTimeSeries(&ts, t);
	printf("X %" NUM_FMT" Y %" NUM_FMT" diff %" NUM_FMT"\n", x, y, x - y);

	t = C(40.0);
	x = InterpolateTimeSeries(interpData, t);
	y = InterpolateTimeSeries(&ts, t);
	printf("X %" NUM_FMT" Y %" NUM_FMT" diff %" NUM_FMT"\n", x, y, x - y);

	t = C(45.0);
	x = InterpolateTimeSeries(interpData, t);
	y = InterpolateTimeSeries(&ts, t);
	printf("X %" NUM_FMT" Y %" NUM_FMT" diff %" NUM_FMT"\n", x, y, x - y);

	t = C(2.0);
	x = InterpolateTimeSeries(interpData, t);
	y = InterpolateTimeSeries(&ts, t);
	printf("X %" NUM_FMT" Y %" NUM_FMT" diff %" NUM_FMT"\n", x, y, x - y);

}

void TestBalance()
{
	int row_count = 100;
	int block_count = 32;
	WetDryRowBound* wet_dry_bounds;
	AllocateWetDryRowBound(row_count, block_count, wet_dry_bounds);



	SGC2_UpdateLoadBalance(100, -1, NULL, wet_dry_bounds);




}

void RunTests()
{






}

//---------------------------------------------------------------------------
// simple linear interpolation from a list of values provided as a boundary condition
NUMERIC_TYPE InterpolateTimeSeries(const NUMERIC_TYPE *varlist, NUMERIC_TYPE t)
{
	int i = 0;
	NUMERIC_TYPE dt, a1, a2;
	NUMERIC_TYPE res = 0;

	// for values less than the start of the array - use 1st value
	if (t < varlist[1]) return(varlist[0]);

	while (varlist[(i + 1) * 2 + 1] > C(-0.9))
	{
		if (varlist[i * 2 + 1] <= t && varlist[i * 2 + 3] > t)
		{
			dt = varlist[i * 2 + 3] - varlist[i * 2 + 1];
			a1 = (t - varlist[i * 2 + 1]) / dt;
			a2 = C(1.0) - a1;

			res = a1*varlist[i * 2 + 2] + a2*varlist[i * 2 + 0];
		}

		i++;
	}

	// for values greater than the end of the array - use last value
	if (t >= varlist[i * 2 + 1])
		res = varlist[i * 2];

	return(res);
}
//-----------------------------------------------------------------------------------

#endif