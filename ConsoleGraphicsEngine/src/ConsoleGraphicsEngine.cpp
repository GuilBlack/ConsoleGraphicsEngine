//#include "Engine3D.h"
#include "vendor/olcConsoleGameEngine.h"
#include <fstream>
#include <strstream>
#include <algorithm>

struct vec3d;

/// <summary>
/// A 4x4 matrix that contains floats.
/// It does some operations with vec3d.
/// </summary>
struct mat4x4f
{
    float mat[4][4] = { 0 };

	#pragma region mat4x4fUtilityMatrices

	static mat4x4f CreateRotationX(float fAngleRad)
	{
		mat4x4f matrix;
		matrix.mat[0][0] = 1.0f;
		matrix.mat[1][1] = cosf(fAngleRad);
		matrix.mat[1][2] = sinf(fAngleRad);
		matrix.mat[2][1] = -sinf(fAngleRad);
		matrix.mat[2][2] = cosf(fAngleRad);
		matrix.mat[3][3] = 1.0f;
		return matrix;
	}

	static mat4x4f CreateRotationY(float fAngleRad)
	{
		mat4x4f matrix;
		matrix.mat[0][0] = cosf(fAngleRad);
		matrix.mat[0][2] = sinf(fAngleRad);
		matrix.mat[2][0] = -sinf(fAngleRad);
		matrix.mat[1][1] = 1.0f;
		matrix.mat[2][2] = cosf(fAngleRad);
		matrix.mat[3][3] = 1.0f;
		return matrix;
	}

	static mat4x4f CreateRotationZ(float fAngleRad)
	{
		mat4x4f matrix;
		matrix.mat[0][0] = cosf(fAngleRad);
		matrix.mat[0][1] = sinf(fAngleRad);
		matrix.mat[1][0] = -sinf(fAngleRad);
		matrix.mat[1][1] = cosf(fAngleRad);
		matrix.mat[2][2] = 1.0f;
		matrix.mat[3][3] = 1.0f;
		return matrix;
	}

	static mat4x4f CreateIdentity()
	{
		mat4x4f matrix;
		matrix.mat[0][0] = 1.0f;
		matrix.mat[1][1] = 1.0f;
		matrix.mat[2][2] = 1.0f;
		matrix.mat[3][3] = 1.0f;
		return matrix;
	}

	static mat4x4f CreateTranslation(float x, float y, float z)
	{
		mat4x4f matrix;
		matrix.mat[0][0] = 1.0f;
		matrix.mat[1][1] = 1.0f;
		matrix.mat[2][2] = 1.0f;
		matrix.mat[3][3] = 1.0f;
		matrix.mat[3][0] = x;
		matrix.mat[3][1] = y;
		matrix.mat[3][2] = z;
		return matrix;
	}

	static mat4x4f CreateProjection(float fFovDegrees, float fAspectRatio, float fNear, float fFar)
	{
		float fFovRad = 1.0f / tanf(fFovDegrees * 0.5f / 180.0f * 3.14159f);
		mat4x4f matrix;
		matrix.mat[0][0] = fAspectRatio * fFovRad;
		matrix.mat[1][1] = fFovRad;
		matrix.mat[2][2] = fFar / (fFar - fNear);
		matrix.mat[3][2] = (-fFar * fNear) / (fFar - fNear);
		matrix.mat[2][3] = 1.0f;
		matrix.mat[3][3] = 0.0f;
		return matrix;
	}

	static mat4x4f CreatePointAt(vec3d& pos, vec3d& target, vec3d& up);

	static mat4x4f QuickInverse(mat4x4f& mat) //only works for rotation/translation matrices that are originally horizontal
	{
		mat4x4f matrix;

		matrix.mat[0][0] = mat.mat[0][0]; matrix.mat[0][1] = mat.mat[1][0]; matrix.mat[0][2] = mat.mat[2][0]; matrix.mat[0][3] = 0.0f;
		matrix.mat[1][0] = mat.mat[0][1]; matrix.mat[1][1] = mat.mat[1][1]; matrix.mat[1][2] = mat.mat[2][1]; matrix.mat[1][3] = 0.0f;
		matrix.mat[2][0] = mat.mat[0][2]; matrix.mat[2][1] = mat.mat[1][2]; matrix.mat[2][2] = mat.mat[2][2]; matrix.mat[2][3] = 0.0f;
		matrix.mat[3][0] = -(mat.mat[3][0] * matrix.mat[0][0] + mat.mat[3][1] * matrix.mat[1][0] + mat.mat[3][2] * matrix.mat[2][0]);
		matrix.mat[3][1] = -(mat.mat[3][0] * matrix.mat[0][1] + mat.mat[3][1] * matrix.mat[1][1] + mat.mat[3][2] * matrix.mat[2][1]);
		matrix.mat[3][2] = -(mat.mat[3][0] * matrix.mat[0][2] + mat.mat[3][1] * matrix.mat[1][2] + mat.mat[3][2] * matrix.mat[2][2]);
		matrix.mat[3][3] = 1.0f;

		return matrix;
	}

	#pragma endregion


	#pragma region mat4x4fOperators

	mat4x4f& operator*=(const mat4x4f& rhs)
	{
		for (int c = 0; c < 4; c++)
			for (int r = 0; r < 4; r++)
				this->mat[r][c] = this->mat[r][0] * rhs.mat[0][c] + this->mat[r][1] * rhs.mat[1][c] + this->mat[r][2] * rhs.mat[2][c] +
				this->mat[r][3] * rhs.mat[3][c];

		return *this;
	}

	mat4x4f operator*(const mat4x4f& rhs) const
	{
		mat4x4f matrix;
		for (int c = 0; c < 4; c++)
			for (int r = 0; r < 4; r++)
				matrix.mat[r][c] = this->mat[r][0] * rhs.mat[0][c] + this->mat[r][1] * rhs.mat[1][c] + this->mat[r][2] * rhs.mat[2][c] +
				this->mat[r][3] * rhs.mat[3][c];

		return matrix;
	}

	#pragma endregion

};

/// <summary>
/// All there is to create a normal vector in 3d space and do some operations on it.
/// </summary>
struct vec3d
{
    float x, y, z, w;

	vec3d()
	{
		x = y = z = 0;
		w = 1;
	}

	vec3d(float a, float b, float c)
	{
		x = a; y = b; z = c; w = 1;
	}

	float GetLength() const
	{
		return sqrtf(x*x + y*y + z*z);
	}

	vec3d Normalize() const
	{
		float l = GetLength();
		return { x / l, y / l, z / l };
	}

	float DotProduct(vec3d& b) const
	{
		return this->x * b.x + this->y * b.y + this->z * b.z;
	}

	vec3d CrossProduct(vec3d& b) const
	{
		vec3d v;
		v.x = this->y * b.z - this->z * b.y;
		v.y = this->z * b.x - this->x * b.z;
		v.z = this->x * b.y - this->y * b.x;
		return v;
	}

	/// <summary>
	/// See if a line intersects with a plane by creating a plane with a point and its normal.
	/// </summary>
	/// <param name="planeP">A point on a plane.</param>
	/// <param name="planeN">A normal associated with the plane.</param>
	/// <param name="lineStart">Start of the line that we want to test against the plane.</param>
	/// <param name="lineEnd">End of the line that we want to test against the plane.</param>
	/// <returns></returns>
	static vec3d LineIntersectsPlane(vec3d& planeP, vec3d& planeN, vec3d& lineStart, vec3d& lineEnd)
	{
		vec3d line = lineEnd - lineStart;
		planeN = planeN.Normalize();
		float planeD = -planeN.DotProduct(planeP);
		float ad = lineStart.DotProduct(planeN);
		float bd = lineEnd.DotProduct(planeN);
		float t = (-planeD - ad) / (bd - ad);
		vec3d lineToIntersect = line * t;
		return lineStart + lineToIntersect;
	}

	#pragma region vec3dOperators

	vec3d& operator+=(const vec3d& rhs)
	{
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
		return *this;
	}

	vec3d& operator-=(const vec3d& rhs)
	{
		this->x -= rhs.x;
		this->y -= rhs.y;
		this->z -= rhs.z;
		return *this;
	}

	vec3d& operator*=(const float rhs)
	{
		this->x *= rhs;
		this->y *= rhs;
		this->z *= rhs;
		return *this;
	}

	vec3d& operator/=(const float rhs)
	{
		this->x /= rhs;
		this->y /= rhs;
		this->z /= rhs;
		return *this;
	}

	vec3d operator+(const vec3d& rhs) const
	{
		return { this->x + rhs.x, this->y + rhs.y, this->z + rhs.z };
	}

	vec3d operator-(const vec3d& rhs) const
	{
		return { this->x - rhs.x, this->y - rhs.y, this->z - rhs.z };
	}

	vec3d operator*(const float rhs) const
	{
		return { this->x * rhs, this->y * rhs, this->z * rhs };
	}

	vec3d operator/(const float rhs) const
	{
		return { this->x / rhs, this->y / rhs, this->z / rhs };
	}

	vec3d operator*(const mat4x4f& rhs) const
	{
		vec3d r;
		r.x = this->x * rhs.mat[0][0] + this->y * rhs.mat[1][0] + this->z * rhs.mat[2][0] + this->w * rhs.mat[3][0];
		r.y = this->x * rhs.mat[0][1] + this->y * rhs.mat[1][1] + this->z * rhs.mat[2][1] + this->w * rhs.mat[3][1];
		r.z = this->x * rhs.mat[0][2] + this->y * rhs.mat[1][2] + this->z * rhs.mat[2][2] + this->w * rhs.mat[3][2];
		r.w = this->x * rhs.mat[0][3] + this->y * rhs.mat[1][3] + this->z * rhs.mat[2][3] + this->w * rhs.mat[3][3];
		return r;
	}

	#pragma endregion

};


/// <summary>
/// Triangles are composed of 3 vertices that are vec3d positions.
/// sym and col represent the symbol used to render on console + color of the symbol.
/// </summary>
struct triangle
{
    vec3d p[3];

	wchar_t sym;
	short col;
};

/// <summary>
/// Stores a vector of triangles that forms an object.
/// </summary>
struct mesh
{
    std::vector<triangle> triangles;

	bool LoadFromObjectFile(std::string sFileName)
	{
		std::ifstream f(sFileName);
		if (!f.is_open())
			return false;

		std::vector<vec3d> vertices;

		while (!f.eof())
		{
			char line[128];
			f.getline(line, 128);

			std::strstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				vec3d v;
				s >> junk >> v.x >> v.y >> v.z;
				vertices.push_back(v);
			}

			if (line[0] == 'f')
			{
				int f[3];
				s >> junk >> f[0] >> f[1] >> f[2];
				triangles.push_back({ vertices[f[0] - 1], vertices[f[1] - 1], vertices[f[2] - 1] });
			}
		}

		return true;
	}
};

#pragma region OtherDeclarations

//Must do the definition here because we didn't define vec3d before mat4x4f
mat4x4f mat4x4f::CreatePointAt(vec3d& pos, vec3d& target, vec3d& up)
{
	//calculating a new forward direction taking into account that our
	//actual forward is from the pos to the target
	vec3d newForward = (target - pos).Normalize();

	//calculating a new up direction
	vec3d a = newForward * up.DotProduct(newForward); //a is a description of the changes between the new forward and the up
	vec3d newUp = (up - a).Normalize();

	//calculating right vector
	vec3d newRight = newUp.CrossProduct(newForward);

	//Constructing Point At matrix
	mat4x4f pointAt;
	pointAt.mat[0][0] = newRight.x;		pointAt.mat[0][1] = newRight.y;		pointAt.mat[0][2] = newRight.z;		pointAt.mat[0][3] = 0.0f;
	pointAt.mat[1][0] = newUp.x;		pointAt.mat[1][1] = newUp.y;		pointAt.mat[1][2] = newUp.z;		pointAt.mat[1][3] = 0.0f;
	pointAt.mat[2][0] = newForward.x;	pointAt.mat[2][1] = newForward.y;	pointAt.mat[2][2] = newForward.z;	pointAt.mat[2][3] = 0.0f;
	pointAt.mat[3][0] = pos.x;			pointAt.mat[3][1] = pos.y;			pointAt.mat[3][2] = pos.z;			pointAt.mat[3][3] = 1.0f;

	return pointAt;

}

#pragma endregion


class Engine3D : public olcConsoleGameEngine
{
public:
    Engine3D()
    {
        m_sAppName = L"3D Graphics Engine";
    }

private:
    mesh m_MeshCube;
    mat4x4f m_MatProj;

	vec3d m_vCamera = { 0, 0, 0 }; //simplified camera
	vec3d m_vLookDir; //where the camera should be looking

	float m_fYaw;


	float m_fTheta;

	CHAR_INFO GetColor(float lum)
	{
		short bg_col, fg_col;
		wchar_t sym;
		int pixel_bw = (int)(13.0f * lum);
		switch (pixel_bw)
		{
		case 0: bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID; break;

		case 1: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_QUARTER; break;
		case 2: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_HALF; break;
		case 3: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 4: bg_col = BG_BLACK; fg_col = FG_DARK_GREY; sym = PIXEL_SOLID; break;

		case 5: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_QUARTER; break;
		case 6: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_HALF; break;
		case 7: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_THREEQUARTERS; break;
		case 8: bg_col = BG_DARK_GREY; fg_col = FG_GREY; sym = PIXEL_SOLID; break;

		case 9:  bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_QUARTER; break;
		case 10: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_HALF; break;
		case 11: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_THREEQUARTERS; break;
		case 12: bg_col = BG_GREY; fg_col = FG_WHITE; sym = PIXEL_SOLID; break;
		default:
			bg_col = BG_BLACK; fg_col = FG_BLACK; sym = PIXEL_SOLID;
		}

		CHAR_INFO c;
		c.Attributes = bg_col | fg_col;
		c.Char.UnicodeChar = sym;
		return c;
	}

	/// <summary>
	/// This function returns an int that indicates the number of triangles to return in the outTris. May be 0, 1, or 2.
	/// </summary>
	/// <param name="planeP">A point in the plane edge that we want to test the triangle against.</param>
	/// <param name="planeN">The normal to the plane edge.</param>
	/// <param name="inTri">The input triangle.</param>
	/// <param name="outTri1">First output triangle if needed.</param>
	/// <param name="outTri2">Second output triangle if needed.</param>
	/// <returns></returns>
	int Triangle_ClipAgainstPlane(vec3d planeP, vec3d planeN, triangle& inTri, triangle& outTri1, triangle& outTri2)
	{
		//normalize the plane normal
		planeN = planeN.Normalize();

		auto dist = [&](vec3d& p)
		{
			vec3d n = p.Normalize();
			//We want to do the dot product between the vector
			//PlaneP -> point of triangle and the unit vec normal to the
			//plane. Basically the same as doing: (n - planeP).planeN
			return (planeN.DotProduct(n) - planeN.DotProduct(planeP));
		};

		vec3d* insidePoints[3]; int nInsidePointCount = 0;
		vec3d* outsidePoints[3]; int nOutsidePointCount = 0;

		float d0 = dist(inTri.p[0]);
		float d1 = dist(inTri.p[1]);
		float d2 = dist(inTri.p[2]);

		if (d0 >= 0) insidePoints[nInsidePointCount++] = &inTri.p[0];
		else outsidePoints[nOutsidePointCount++] = &inTri.p[0];
		if (d1 >= 0) insidePoints[nInsidePointCount++] = &inTri.p[1];
		else outsidePoints[nOutsidePointCount++] = &inTri.p[1];
		if (d2 >= 0) insidePoints[nInsidePointCount++] = &inTri.p[2];
		else outsidePoints[nOutsidePointCount++] = &inTri.p[2];

		if (nInsidePointCount == 0)
		{
			//all points are outside the plane so we want to clip the whole triangle.
			//it will cease to exist.
			return 0; //no return triangles are valid.
		}

		if (nInsidePointCount == 3)
		{
			//all points are inside the plane so we want to return the whole triangle.
			outTri1 = inTri;
			return 1; //no return triangles are valid.
		}

		if (nInsidePointCount == 1 && nOutsidePointCount == 2)
		{
			//the triangle should be clipped. Since only 1 point is inside,
			//the triangle will become smaller. No need for subdivision.

			//outTri1.col = inTri.col;
			outTri1.col = FG_GREEN;
			outTri1.sym = inTri.sym;

			outTri1.p[0] = *insidePoints[0];

			outTri1.p[1] = vec3d::LineIntersectsPlane(planeP, planeN, *insidePoints[0], *outsidePoints[0]);
			outTri1.p[2] = vec3d::LineIntersectsPlane(planeP, planeN, *insidePoints[0], *outsidePoints[1]);

			return 1;
		}

		if (nInsidePointCount == 2 && nOutsidePointCount == 1)
		{
			//the triangle should be clipped. Since 2 points are inside,
			//The triangle will be subdivided into 2 smaller triangles.
			//outTri1.col = inTri.col;
			outTri1.col = FG_BLUE;
			outTri1.sym = inTri.sym;

			//outTri2.col = inTri.col;
			outTri2.col = FG_CYAN;
			outTri2.sym = inTri.sym;

			outTri1.p[0] = *insidePoints[0];
			outTri1.p[1] = *insidePoints[1];
			outTri1.p[2] = vec3d::LineIntersectsPlane(planeP, planeN, *insidePoints[0], *outsidePoints[0]);

			outTri2.p[0] = *insidePoints[1];
			outTri2.p[1] = outTri1.p[2];
			outTri2.p[2] = vec3d::LineIntersectsPlane(planeP, planeN, *insidePoints[1], *outsidePoints[0]);

			return 2;
		}
	}

public:
    bool OnUserCreate() override
	{

		m_MeshCube.LoadFromObjectFile("res/obj/UtahTeapot.obj");

		// Projection Matrix
		m_MatProj = mat4x4f::CreateProjection(90, (float)ScreenHeight() / (float)ScreenWidth(), 0.1f, 1000.0f);

		return true;
	}

    bool OnUserUpdate(float fElapseTime) override
	{
		if (GetKey(VK_UP).bHeld)
			m_vCamera.y -= 8.0f * fElapseTime;

		if (GetKey(VK_DOWN).bHeld)
			m_vCamera.y += 8.0f * fElapseTime;
		
		if (GetKey(VK_RIGHT).bHeld)
			m_vCamera.x += 8.0f * fElapseTime;

		if (GetKey(VK_LEFT).bHeld)
			m_vCamera.x -= 8.0f * fElapseTime;

		vec3d vForward = m_vLookDir * (8.0f * fElapseTime); //creates a velocity vector going in the direction that the cam is facing

		if (GetKey(L'W').bHeld)
			m_vCamera += vForward;

		if (GetKey(L'S').bHeld)
			m_vCamera -= vForward;

		if (GetKey(L'A').bHeld)
			m_fYaw += 2.0f * fElapseTime;

		if (GetKey(L'D').bHeld)
			m_fYaw -= 2.0f * fElapseTime;

		Fill(0, 0, ScreenWidth(), ScreenHeight(), PIXEL_SOLID, FG_BLACK);

		//Set Rotation Matrices
		mat4x4f matRotZ, matRotX;
		//m_fTheta += 1.0f * fElapseTime;

		matRotZ = mat4x4f::CreateRotationZ(3.1415f);
		matRotX = mat4x4f::CreateRotationX(m_fTheta);

		//Set Translation Matrix
		mat4x4f matTrans;
		matTrans = mat4x4f::CreateTranslation(0.0f, 1.0f, 5.0f);

		//Set World Matrix
		mat4x4f matWorld = mat4x4f::CreateIdentity();
		matWorld = matRotZ * matRotX;
		matWorld *= matTrans;

		vec3d vUp = { 0, 1, 0 };
		vec3d vTarget = { 0, 0, 1 }; //takes target vector that looks forward on the zAxis
		mat4x4f matCameraRot = mat4x4f::CreateRotationY(m_fYaw);
		m_vLookDir = vTarget * matCameraRot; //rotates the target vector to create a look dir vector
		vTarget = m_vCamera + m_vLookDir; //offsetting the target be the camera + where it's looking at

		mat4x4f matCamera = mat4x4f::CreatePointAt(m_vCamera, vTarget, vUp);

		mat4x4f matView = mat4x4f::QuickInverse(matCamera);

		std::vector<triangle> trianglesToRaster;

		//Draw Triangles
		for (const auto &tri : m_MeshCube.triangles)
		{
			//rotates the cube
			triangle triProjected, triTransformed, triViewed;

			triTransformed.p[0] = tri.p[0] * matWorld;
			triTransformed.p[1] = tri.p[1] * matWorld;
			triTransformed.p[2] = tri.p[2] * matWorld;

			vec3d normal, vecA, vecB; //normal is the normal vector to the plane made with the 2 vectors a and b
			vecA = triTransformed.p[1] - triTransformed.p[0];
			vecB = triTransformed.p[2] - triTransformed.p[0];

			normal = vecA.CrossProduct(vecB); //cross prod between the 2 lines to get the normal
			normal = normal.Normalize(); //normalizing the normal

			//vector between the surface normal and the camera
			vec3d projLine = triTransformed.p[0] - m_vCamera;

			//draw triangle only if the normal is facing the cam so if its z is negative
			if (projLine.DotProduct(normal) < 0)
			{
				//Lighting
				vec3d directionalLight = { 0.0f, 0.0f, -1.0f }; //light facing the the cube
				directionalLight = directionalLight.Normalize();

				float dotProdNormLight = normal.DotProduct(directionalLight);

				CHAR_INFO c = GetColor(dotProdNormLight);
				triTransformed.col = c.Attributes;
				triTransformed.sym = c.Char.UnicodeChar;

				// Converts world space to the camera's view
				triViewed.p[0] = triTransformed.p[0] * matView;
				triViewed.p[1] = triTransformed.p[1] * matView;
				triViewed.p[2] = triTransformed.p[2] * matView;

				//Clip Viewed Triangle against near plane.
				//Could form at most 2 triangles.
				int nClippedTriangles = 0;
				triangle clipped[2];
				nClippedTriangles = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.1f }, { 0.0f, 0.0f, 1.0f }, triViewed, clipped[0], clipped[1]);

				for (int n = 0; n < nClippedTriangles; n++)
				{
					// Projects the triangles from 3d to 2d
					triProjected.p[0] = clipped[n].p[0] * m_MatProj;
					triProjected.p[1] = clipped[n].p[1] * m_MatProj;
					triProjected.p[2] = clipped[n].p[2] * m_MatProj;

					triProjected.col = triTransformed.col;
					triProjected.sym = triTransformed.sym;

					// Scale into view, we normalize the coordinates into a 2d-kind-of space by dividing the vector by the depth z.
					triProjected.p[0] = triProjected.p[0] / triProjected.p[0].w;
					triProjected.p[1] = triProjected.p[1] / triProjected.p[1].w;
					triProjected.p[2] = triProjected.p[2] / triProjected.p[2].w;

					//Scale into view

					//will make range -1 / 1 to 0 / 2
					vec3d vOffsetView = { 1, 1, 0 };
					triProjected.p[0] += vOffsetView;
					triProjected.p[1] += vOffsetView;
					triProjected.p[2] += vOffsetView;

					//rescale it between 0 and 1 then adjusting it to screen sizes
					triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
					triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
					triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
					triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
					triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
					triProjected.p[2].y *= 0.5f * (float)ScreenHeight();

					trianglesToRaster.push_back(triProjected);
				}
			}

		}

		//based on the painter algorithm to draw the furthest triangles first
		std::sort(trianglesToRaster.begin(), trianglesToRaster.end(), 
			[](triangle &t1, triangle &t2)
			{
				//taking the z of the mid point of the faces and comparing them
				float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
				float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
				return z1 > z2;
			});

		for (auto& triToRaster : trianglesToRaster)
		{
			triangle clipped[2];
			std::list<triangle> listTriangles;
			listTriangles.push_back(triToRaster);
			int nNewTriangles = 1;

			for (int p = 0; p < 4; p++)
			{
				int nTrisToAdd = 0;

				//will do so until all the known triangles
				//in this itteration where we test against
				//a spefic plane are clipped.
				while (nNewTriangles > 0)
				{
					//takes the triangle in front of the queue
					triangle triToTest = listTriangles.front();
					listTriangles.pop_front();
					--nNewTriangles;

					//tests the triangles against each planes.
					//if one or more triangle(s) is/are returned,
					//they will be added to the list and tested
					//against the other planes.
					switch (p)
					{
					case 0: nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, triToTest, clipped[0], clipped[1]); break;
					case 1: nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, (float)ScreenHeight() - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, triToTest, clipped[0], clipped[1]); break;
					case 2: nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, triToTest, clipped[0], clipped[1]); break;
					case 3: nTrisToAdd = Triangle_ClipAgainstPlane({ (float)ScreenWidth() - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, triToTest, clipped[0], clipped[1]); break;
					}

					//add the newly formed triangles if any.
					for (int t = 0; t < nTrisToAdd; t++)
						listTriangles.push_back(clipped[t]);
				}
				nNewTriangles = listTriangles.size();
			}

			for (auto& t : listTriangles)
			{
				FillTriangle((int)t.p[0].x, (int)t.p[0].y, (int)t.p[1].x, (int)t.p[1].y, (int)t.p[2].x, (int)t.p[2].y, t.sym, t.col);

				//DrawTriangle((int)t.p[0].x, (int)t.p[0].y, (int)t.p[1].x, (int)t.p[1].y, (int)t.p[2].x, (int)t.p[2].y, PIXEL_SOLID, FG_RED);

			}

		}

		return true;
	}
};

int main()
{
    Engine3D engine;
	if (engine.ConstructConsole(256, 156, 4, 4))
		engine.Start();
	else
		std::cout << "[ERROR]: Console didn't start" << std::endl;

	return 0;
}
