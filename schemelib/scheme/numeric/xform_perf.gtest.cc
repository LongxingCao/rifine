#include <gtest/gtest.h>

#include "scheme/numeric/rand_xform.hh"

#include <Eigen/Geometry>

#include <random>
#include "scheme/util/Timer.hh"

namespace scheme { namespace numeric { namespace pref_test {

using std::cout;
using std::endl;

typedef Eigen::Transform<double,3,Eigen::AffineCompact> Xform;
// typedef Eigen::Affine3d Xform;


// TEST( XformMap, basic_test ){
template< class Xform >
void test_xform_perf(){
	int NSAMP = 100*1000;
	#ifdef SCHEME_BENCHMARK
	NSAMP = 10*1000*1000;
	#endif
	std::mt19937 rng;
	Xform x,sum = Xform::Identity();
	rand_xform(rng,x);

	util::Timer<> t;
	for(int i = 0; i < NSAMP; ++i){
		sum = sum * x;
	}
	double time = t.elapsed_nano();
	printf( "runtime %7.3fns nonsense: %7.3f \n", time/NSAMP, sum.translation()[0] );

}

TEST( xform_perf, preformance ){
	cout << "AffineCompact d "; test_xform_perf< Eigen::Transform<double,3,Eigen::AffineCompact> >();
	cout << "Affine        d "; test_xform_perf< Eigen::Affine3d  >();
	cout << "AffineCompact f "; test_xform_perf< Eigen::Transform<float ,3,Eigen::AffineCompact> >();
	cout << "Affine        f "; test_xform_perf< Eigen::Affine3f  >();
	// cout << "XformHash_bt24_Cubic_Zorder"; test_xform_perf< XformHash_bt24_Cubic_Zorder >();
}

}}}
