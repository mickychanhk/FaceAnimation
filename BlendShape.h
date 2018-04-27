#pragma once

#include "common.h"

class BlendShape
{
public:
	struct vert_t {
		float x, y, z;
	};
	typedef vector<vert_t> shape_t;

	BlendShape(void);
	~BlendShape(void);

	const shape_t& expression(int idx) const {
		return exprList[idx];
	}

	int expressionCount() const {
		return nShapes;
	}	

public:
	bool read(const string& filename);

private:
	int nVerts;
	int nShapes;
	int nFaces;

	vector<shape_t> exprList;
};

