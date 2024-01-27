
test:
	clib install --dev
	@$(CC) test.c -std=c99 -I src -I deps -o $@ -lm
	@./$@

.PHONY: test
