#ifndef INCLUDED_chemical_stub_HH
#define INCLUDED_chemical_stub_HH

#include <Eigen/Geometry>

namespace scheme {
namespace chemical {


    template<class Xform,class InputPoint>
    Xform make_stub(InputPoint cen, InputPoint n, InputPoint ca, InputPoint c){
        typedef Eigen::Matrix< typename Xform::Scalar, 3, 1 > XPoint;
        XPoint _n, _ca, _c, _cen;
        for( int i = 0; i < 3; ++i){
            _n  [i] = n  [i];
            _ca [i] = ca [i];
            _c  [i] = c  [i];
            _cen[i] = cen[i];
        }
        Xform out;
        XPoint e1( _n - _ca );
        e1.normalize();
        XPoint e3( e1.cross(_c-_ca) );
        e3.normalize();
        XPoint e2( e3.cross(e1) );
        out.linear().col(0) = e1;
        out.linear().col(1) = e2;
        out.linear().col(2) = e3;
        out.translation() = _cen;
        return out;
    }

    template<class Xform,class InputPoint>
    Xform make_stub(InputPoint n, InputPoint ca, InputPoint c){
        return make_stub<Xform>( ca, n, ca, c );
    }

}
}

#endif
