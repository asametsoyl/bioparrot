import { defineConfig } from 'vitest/config';
import { WxtVitest } from 'wxt/testing/vitest-plugin';

export default defineConfig({
  plugins: [WxtVitest()],
  test: {
    environment: 'happy-dom',
    include: process.env.BP_PDF ? ['tests/debug/**/*.test.ts'] : ['tests/unit/**/*.test.ts'],
    exclude: ['**/node_modules/**', '**/._*'],
  },
});
