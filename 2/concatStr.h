#ifndef CONCATSTR_H
#define CONCATSTR_H

#include <stdlib.h>
#include <string.h>

inline char* concatStr(const char* s1, const char* s2) {
	char* result = (char*)malloc(sizeof(char) * (strlen(s1) + strlen(s2) + 1));

	strcpy(result, s1);
	strcat(result, s2);

	return result;
}

#endif