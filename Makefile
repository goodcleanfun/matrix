
test:
	clib install --dev --concurrency 1
	@$(CC) test.c deps/file_utils/file_utils.c -std=c99 -I src -I deps -o $@ -lm
	@./$@

.PHONY: test
