/***************************************************************************
 *   Copyright (C) 2004 by Olivier Galizzi                                 *
 *   olivier.galizzi@imag.fr                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Mesh2D.hpp"
#include "OpenGLWrapper.hpp"
#include <iostream>


Mesh2D::Mesh2D () : BodyGeometricalModel()
{		
	createIndex();
}

Mesh2D::~Mesh2D ()
{

}

void Mesh2D::postProcessAttributes(bool deserializing)
{
	BodyGeometricalModel::postProcessAttributes(deserializing);
	
	if(deserializing)
	{
	
		fNormals.resize(faces.size());
		vNormals.resize(vertices.size());
		triPerVertices.resize(vertices.size());
		for(unsigned int i=0;i<faces.size();i++)
		{
			triPerVertices[faces[i][0]].push_back(i);
			triPerVertices[faces[i][1]].push_back(i);
			triPerVertices[faces[i][2]].push_back(i);
		}
	}
}

void Mesh2D::registerAttributes()
{
	BodyGeometricalModel::registerAttributes();
	REGISTER_ATTRIBUTE(vertices);
	REGISTER_ATTRIBUTE(edges);
	REGISTER_ATTRIBUTE(faces);
	REGISTER_ATTRIBUTE(width);
	REGISTER_ATTRIBUTE(height);
}



