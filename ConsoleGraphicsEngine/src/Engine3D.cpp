//#include "Engine3D.h"
//
//Engine3D::Engine3D()
//{
//	m_sAppName = L"3D Graphics Engine";
//}
//
//void Engine3D::MultiplyMatrixVector(const vec3d& vecIn, vec3d vecOut, const mat4x4f& mat)
//{
//	vecOut.x = vecIn.x * mat.mat[0][0] + vecIn.y * mat.mat[1][0] + vecIn.z * mat.mat[2][0] + mat.mat[3][0];
//	vecOut.y = vecIn.x * mat.mat[0][1] + vecIn.y * mat.mat[1][1] + vecIn.z * mat.mat[2][1] + mat.mat[3][1];
//	vecOut.z = vecIn.x * mat.mat[0][2] + vecIn.y * mat.mat[1][2] + vecIn.z * mat.mat[2][2] + mat.mat[3][2];
//	float w = vecIn.x * mat.mat[0][2] + vecIn.y * mat.mat[1][2] + vecIn.z * mat.mat[2][2] + mat.mat[3][2];
//
//	if (w != 0.0f)
//	{
//		vecOut.x /= w; vecOut.y /= w; vecOut.z /= w;
//	}
//}
//
//bool Engine3D::OnUserCreate()
//{
//	m_MeshCube.triangles = {
//		// SOUTH
//		{ 0.0f, 0.0f, 0.0f,    0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 0.0f },
//		{ 0.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 0.0f, 0.0f },
//
//		// EAST                                                      
//		{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f },
//		{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 0.0f, 1.0f },
//
//		// NORTH                                                     
//		{ 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f },
//		{ 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f },
//
//		// WEST                                                      
//		{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f },
//		{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f,    0.0f, 0.0f, 0.0f },
//
//		// TOP                                                       
//		{ 0.0f, 1.0f, 0.0f,    0.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f },
//		{ 0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 0.0f },
//
//		// BOTTOM                                                    
//		{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f },
//		{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f,    1.0f, 0.0f, 0.0f },
//	};
//
//	// Projection Matrix
//	float fNear = 0.1f;
//	float fFar = 1000.0f;
//	float fFov = 75.0f;
//	float fAspectRatio = (float)ScreenHeight() / (float)ScreenWidth();
//	float fFovRad = 1.0f / tanf(fFov * 0.5f / 180.0f * 3.14159265f); //in radian and not in degrees
//
//	m_MatProj.mat[0][0] = fAspectRatio * fFovRad;
//	m_MatProj.mat[1][1] = fFovRad;
//	m_MatProj.mat[2][2] = fFar / (fFar - fNear);
//	m_MatProj.mat[3][3] = (-fNear * fFar) / (fFar - fNear);
//	m_MatProj.mat[2][3] = 1.0f;
//	m_MatProj.mat[3][3] = 0.0f;
//
//	return true;
//}
//
//bool Engine3D::OnUserUpdate(float fElapseTime)
//{
//	Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);
//
//	//Draw Triangles
//	for (auto tri : m_MeshCube.triangles)
//	{
//		triangle triProjected;
//		MultiplyMatrixVector(tri.p[0], triProjected.p[0], m_MatProj);
//		MultiplyMatrixVector(tri.p[1], triProjected.p[1], m_MatProj);
//		MultiplyMatrixVector(tri.p[2], triProjected.p[2], m_MatProj);
//
//		DrawTriangle(triProjected.p[0].x, triProjected.p[0].y,
//			triProjected.p[1].x, triProjected.p[1].y,
//			triProjected.p[2].x, triProjected.p[2].y,
//			PIXEL_SOLID, FG_WHITE);
//	}
//
//	return true;
//}
