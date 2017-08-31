#ifndef __Point_h
#define __Point_h

/**
	class name: Point
	Nested Point class
	Represents one point in dataset.*/
	class Point
	{
	public:
		double	DCoord[3];
		double  DNormal[3];
		bool	interpolated;

		int		Id;
		int		TendonConnectionID; // for border Fiber PointIDs, if BBoundary, this ID corresponds to tendon
		bool	BMuscle;	// indicates, if point is part of muscle or tendon
		bool    BBoundary;
		static int	count;

		// for interpolation purposes:
		double	a[3], b[3], c[3], d[3]; // (3 axes)
		double	distToNext;
		double  t;

		double thickness;

		/** constructor */
		Point(double *pCoord);
		/** destructor */
		~Point();

		void ComputeInterpolatedCurveCoord(double t, double * result, bool substract);
		// Assigns normal to this Point
		void SetNormal(double* normal);
	};

#endif