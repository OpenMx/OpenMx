#ifndef _FREEBSDMERGE_H_
#define _FREEBSDMERGE_H_

typedef int (*mergesort_cmp_t)(const void *, const void *, void *userdata);

int freebsd_mergesort(void *base, size_t nmemb, size_t size, mergesort_cmp_t cmp, void *userdata);


#endif // _FREEBSDMERGE_H_
