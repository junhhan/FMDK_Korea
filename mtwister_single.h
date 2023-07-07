/* 
   A C-program for MT19937, with initialization improved 2002/2/10.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   This is a faster version by taking Shawn Cokus's optimization,
   Matthe Bellew's simplification, Isaku Wada's real version.

   Before using, initialize the state by using init_genrand(seed) 
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/CODES/MTARCOK/mt19937ar-cok.c 
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/
#define MTn 624
#define MTm 397
#define MATRIX_A 0x9908b0dfUL
#define UMASK 0x80000000UL
#define LMASK 0x7fffffffUL
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))

static unsigned long MTstate[MTn]; // an array for the state vector
static int MTleft= 1;
static int MTinitf= 0;
static unsigned long *MTnext;

void _MTinit(unsigned long seed) {
	int i;
	MTstate[0]= seed & 0xffffffffUL;
	for (i= 1; i < MTn; ++i) {
		MTstate[i] = (1812433253UL * (MTstate[i-1] ^ (MTstate[i-1] >> 30)) + i); 
		MTstate[i] &= 0xffffffffUL;
	}
	MTleft= 1; MTinitf= 1;
}

static void _MTnext(void) {
	unsigned long *p= MTstate;

	if (MTinitf == 0) {
		_MTinit(5489UL);
	}
	
	MTleft= MTn;
    MTnext= MTstate;
	
	int i;
	for (i= MTn - MTm + 1; --i; ++p) {
		*p= p[MTm] ^ TWIST(p[0], p[1]);
	} 
	
	for (i= MTm; --i; ++p) {
		*p = p[MTm - MTn] ^ TWIST(p[0], p[1]);
	} 
	
	*p= p[MTm - MTn] ^ TWIST(p[0], MTstate[0]);
}

double _Rand(void) {
    unsigned long y;

    MTleft -=1;
	if (MTleft == 0) {_MTnext();}
    y= *MTnext++;

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return ((double)y * (1.0/4294967295.0));
}

unsigned long _Rand_ulong(void) {
	unsigned long y;

    MTleft -=1;
	if (MTleft == 0) {_MTnext();}
    y= *MTnext++;

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}
