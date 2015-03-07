/* 
 * BOUNDING SPHERE v1.0
 * T.J. Moran (http://tjm.io)
 *
 * This software is dedicated to the public domain.
 */

#include <math.h>
#include <stdlib.h>
#include <float.h>

/* bounding_sphere()
 *
 * Computes a bounding sphere for a list of points. The result is approximate in
 * the sense that it may not be the smallest possible bounding sphere, but not 
 * in the sense that it may not bound the points; it always will.
 * 
 * The procedure used is the iterative variant of [1] described in [2].
 * 
 * The input is an array of vertices, each beginning with a vector of 'ndim'
 * floats. The 'bytes_per_vertex' parameter gives the size of each vertex. 
 * Vertices may be larger than the vectors they contain, allowing typical 
 * polymesh vertex buffers to be passed directly to the function.
 * 
 * For portable behaviour the caller should ensure that 'bytes_per_vertex' is
 * a multiple of sizeof(float) and that the input array is likewise aligned. The
 * code will generate unaligned loads if these conditions are not met, which do
 * not work on all platforms (e.g. many ARM CPUs).
 * 
 * [1] Ritter, J. An efficient bounding sphere. In Andrew S. Glassner, editor, 
 *     Graphics Gems. Academic Press, Boston, MA, 1990.
 * [2] Ericson, C. Real-Time Collision Detection. Morgan Kaufmann, 2004.
 * 
 * RETURN VALUE
 * 
 *     The radius of the bounding sphere.
 *
 * PARAMETERS
 *
 *     out_center       - Location to write the center point of the sphere.
 *     vertices         - An array of vertices, each beginning with an N-d vector.
 *     num_vertices     - Number of entries in the vertex array.
 *     bytes_per_vertex - A stride in bytes used to access the vertex array.
 *     ndim             - Number of components in the vector to read from each 
 *                        vertex.
 *     indices          - Optional array of 'num_vertices' indices into the 
 *                        input array. The indices are vertex numbers, not byte 
 *                        offsets. Useful for processing a subset of a larger
 *                        vertex array.
 *     num_iterations   - Number of times to repeat the expansion phase in an 
 *                        attempt to improve the solution. For typical inputs 
 *                        there is usually little benefit in more than 4 or 8 
 *                        iterations.
 */
float bounding_sphere(
	float *out_center,
	const void *vertices, 
	unsigned num_vertices, 
	unsigned bytes_per_vertex,
	unsigned ndim,
	const unsigned *indices,
	unsigned num_iterations);

#if defined(BSPHERE_IMPLEMENTATION)

float bounding_sphere(
	float *out_center,
	const void *vertices, 
	unsigned num_vertices, 
	unsigned bytes_per_vertex,
	unsigned ndim,
	const unsigned *indices,
	unsigned num_iterations)
{
	static const float SHRINKING_FACTOR = 0.95f;

	/* The result is undefined if there are no vertices. */
	if (num_vertices == 0 || bytes_per_vertex < ndim * sizeof(float)) {
		if (out_center != NULL) {
			for (unsigned i = 0; i < ndim; ++i)
				out_center[i] = 0.0f;
		}
		return 0.0f;
	}

	/* Allocate a buffer for the permutation array and the N-dimensional 
	 * temporary vectors we need. */
	unsigned pvsize = ndim * sizeof(float *);
	unsigned vsize  = ndim * sizeof(float);
	unsigned buffer_size = num_vertices * sizeof(float *); /* Permutation. */
	buffer_size += 2 * pvsize; /* Extreme point arrays. */
	buffer_size += 3 * vsize; /* Center vectors and radius vector. */
	char *buf = (char *)malloc(buffer_size);
	const float **min = (const float **)(buf + 0 * pvsize);
	const float **max = (const float **)(buf + 1 * pvsize);
	float *center      = (float *)(buf + 2 * pvsize + 0 * vsize);
	float *best_center = (float *)(buf + 2 * pvsize + 1 * vsize);
	float *rv          = (float *)(buf + 2 * pvsize + 2 * vsize);
	const float **permutation = (const float **)(buf + 2 * pvsize + 3 * vsize);

	/* Populate the permutation array with pointers to the input vectors. */
	const char *vertex_bytes = (const char *)vertices;
	if (indices != NULL) {
		for (unsigned i = 0; i < num_vertices; ++i)
			permutation[i] = (float *)(vertex_bytes + indices[i] * bytes_per_vertex);
	} else {
		for (unsigned i = 0; i < num_vertices; ++i)
			permutation[i] = (float *)(vertex_bytes + i * bytes_per_vertex);
	}

	/* Find the pair of extreme points (min, max) on each axis. */
	for (unsigned i = 0; i < ndim; ++i)
		min[i] = max[i] = permutation[0];
	for (unsigned i = 1; i < num_vertices; ++i) {
		const float *vertex = permutation[i];
		for (unsigned j = 0; j < ndim; ++j) {
			if (vertex[j] < min[j][j])
				min[j] = vertex;
			if (vertex[j] > max[j][j]) 
				max[j] = vertex;
		}
	}

	/* The support pair is the furthest-apart pair of extreme points. */
	unsigned support = 0;
	float support_dsq = 0.0f;
	for (unsigned i = 0; i < ndim; ++i) {
		float dsq = 0.0f;
		for (unsigned j = 0; j < ndim; ++j) {
			float d = max[i][j] - min[i][j];
			dsq += d * d;
		}
		if (dsq > support_dsq) {
			support = i;
			support_dsq = dsq;
		}
	}

	/* Make a sphere enclosing the support pair. */
	float radius = 0.0f, best_radius = FLT_MAX;
	for (unsigned i = 0; i < ndim; ++i) {
		center[i] = 0.5f * (min[support][i] + max[support][i]);
		float d = center[i] - min[support][i];
		radius += d * d;
	}
	float rr = sqrtf(radius);
	if (rr != 0.0f)
		radius /= rr;

	/* Expand the sphere to enclose all points still outside it, using a 
	 * different vertex order in each iteration in the hope of improving the
	 * result. */
	do {
		float radius_squared = radius * radius;
		for (unsigned i = 0; i < num_vertices; ++i) {
			/* Swap a random, later element of the permutation array with this 
			 * one. */
			unsigned j = rand() % (num_vertices - i);
			const float *other = permutation[i + j];
			permutation[i + j] = permutation[i];
			permutation[i] = other;

			/* Compute the vector from the center to the current point and its
			 * squared length. */
			float dsq = 0.0f;
			for (unsigned j = 0; j < ndim; ++j) {
				rv[j] = other[j] - center[j];
				dsq += rv[j] * rv[j];
			}

			/* Expand the sphere to enclose the point. */
			if (dsq > radius_squared) {
				float distance = sqrtf(dsq);
				float new_radius = (radius + distance) * 0.5f;
				float radius_delta = new_radius - radius;
				radius = new_radius;
				radius_squared = new_radius * new_radius;
				float s = radius_delta / distance;
				for (unsigned j = 0; j < ndim; ++j)
					center[j] += rv[j] * s;
			}
		}

		/* Is the result better than the current best? */
		if (radius < best_radius) {
			for (unsigned i = 0; i < ndim; ++i)
				best_center[i] = center[i];
			best_radius = radius;
		}

		/* Shrink the sphere and try again with a different vertex order. */
		radius *= SHRINKING_FACTOR;
	} while (num_iterations--);
	
	if (out_center != NULL) {
		for (unsigned i = 0; i < ndim; ++i)
			out_center[i] = best_center[i];
	}

	free(buf);

	return best_radius;
}

#endif /* defined(BSPHERE_IMPLEMENTATION) */
