#ifndef INCLUDED_kinematics_Director_meta.HH
#define INCLUDED_kinematics_Director_meta.HH

#include <scheme/util/meta/util.hh>
#include <scheme/types.hh>


namespace scheme {
namespace kinematics {
template<
    class _Position,
    class _BigIndex = uint64_t,
    class _Index = uint64_t
>
struct Director;
}
}

using ::scheme::make_shared;
using ::scheme::shared_ptr;

template <class __Director>
using _DirectorPosition = typename ::scheme::util::meta::remove_pointer<__Director>::type::Position;
template <class __Director>
using _DirectorIndex = typename ::scheme::util::meta::remove_pointer<__Director>::type::Index;
template <class __Director>
using _DirectorBigIndex = typename ::scheme::util::meta::remove_pointer<__Director>::type::BigIndex;
template <class __Director>
using _DirectorScaffoldIndex = typename ::scheme::util::meta::remove_pointer<__Director>::type::ScaffoldIndex;

template <class __Director>
using _DirectorBase = shared_ptr< ::scheme::kinematics::Director<_DirectorPosition<__Director>,
    _DirectorBigIndex<__Director>,
    _DirectorIndex<__Director>
     >>;

template <class __Director>
using _ScaffoldIndexFromDirector = typename _DirectorBigIndex<__Director>::type::second_type;

#endif
