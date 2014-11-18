#pragma once

namespace Random {
	/**
	 * Pseudo Random Number Generator
	 * @return a Pseudo-Random Number from 0 to 2^27 -1
	 */
	long PRNG();

	void seed(long a_seed);

	/**
	 * @return relatively random bit sequence (on a fast CPU). Use sparingly!
	 * Expected runtime is about (a_numBits*2)+1 milliseconds (a little less),
	 * during which 100% of CPU is used
	 */
	long TRNG(int a_numBits);

	float PRNGf();

	float PRNGf(float min, float max);
};
