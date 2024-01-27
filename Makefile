
test:
	clib install --dev
	@$(CC) test.c deps/file_utils/file_utils.c -std=c99 -I src -I deps -o $@ -lm
	@./$@

.PHONY: test
