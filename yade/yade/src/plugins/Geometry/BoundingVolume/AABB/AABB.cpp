#include "AABB.hpp"
#include "OpenGLWrapper.hpp"

AABB::AABB () : BodyBoundingVolume()
{
	createIndex();
}

void AABB::postProcessAttributes(bool /*deserializing*/)
{
// 	if(deserializing)
// 	{
// 		initialHalfSize = halfSize;
// 		initialCenter = center;
// 		min = center-halfSize;
// 		max = center+halfSize;
// 	}
}

void AABB::registerAttributes()
{
	BodyBoundingVolume::registerAttributes();
	//REGISTER_ATTRIBUTE(center);
	//REGISTER_ATTRIBUTE(halfSize);
}




//FIXME : remove all updates and overlap function from BodyBoundingVolume
// void AABB::update(Se3r& se3)
// {
// 	//center = se3.translation;
// 	Vector3r newHalf,newCenter;
// 
//  	center = se3.translation;// + se3.rotation * initialCenter;
// 
// 	Matrix3r r;
// 	se3.rotation.toRotationMatrix(r);
// 	halfSize = Vector3r(0,0,0);
// 	for( int i=0; i<3; ++i )
// 		for( int j=0; j<3; ++j )
// 			halfSize[i] += fabs( r[i][j] * initialHalfSize[j] );
// 
// 	min = center-halfSize;
// 	max = center+halfSize;
// 
// }
// 
// bool AABB::overlap(const AABB& aabb)
// {
// 	Vector3r max2 = aabb.max;
// 	Vector3r min2 = aabb.min;
// 
// 	return !(max[0]<min2[0] || max2[0]<min[0] || max[1]<min2[1] || max2[1]<min[1] || max[2]<min2[2] || max2[2]<min[2]);
// }
