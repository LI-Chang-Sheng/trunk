#ifndef __AABB_H__
#define __AABB_H__

#include "BodyBoundingVolume.hpp"

class AABB : public BodyBoundingVolume
{
	public : Vector3r halfSize;
	public : Vector3r center;

	private : Vector3r initialHalfSize;
	private : Vector3r initialCenter;


	// construction
	public : AABB ();

///////////////////////////////////////////////////////////////////////////////////////////////////
/// Serialization										///
///////////////////////////////////////////////////////////////////////////////////////////////////
	
	REGISTER_CLASS_NAME(AABB);	
	protected : virtual void postProcessAttributes(bool deserializing);
	public : void registerAttributes();
	
///////////////////////////////////////////////////////////////////////////////////////////////////
/// Indexable											///
///////////////////////////////////////////////////////////////////////////////////////////////////
	
	REGISTER_CLASS_INDEX(AABB,BodyBoundingVolume);
};

REGISTER_SERIALIZABLE(AABB,false);

#endif // __AABB_H__



