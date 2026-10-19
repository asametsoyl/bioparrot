import { defineConfig } from 'wxt';

// https://wxt.dev/api/config.html
export default defineConfig({
  srcDir: 'src',
  modules: ['@wxt-dev/module-react'],
  imports: false,
  zip: {
    exclude: ['**/._*', '**/.DS_Store'],
    excludeSources: ['**/._*', '**/.DS_Store', '**/node_modules/**', '.output/**', 'test-results/**'],
  },
  manifest: ({ browser }) => ({
    name: '__MSG_extName__',
    short_name: 'BioParrot',
    description: '__MSG_extDescription__',
    default_locale: 'tr',
    permissions: ['storage', 'unlimitedStorage', 'contextMenus', 'scripting'],
    host_permissions: ['<all_urls>'],
    action: {
      default_title: '__MSG_extName__',
    },
    // PDF sayfalarından BioParrot PDF çevirmenine yönlendirme için
    web_accessible_resources: [{ resources: ['pdf.html'], matches: ['<all_urls>'] }],
    commands: {
      'toggle-translate': {
        suggested_key: { default: 'Alt+A' },
        description: '__MSG_cmdToggleTranslate__',
      },
      'toggle-display-mode': {
        suggested_key: { default: 'Alt+S' },
        description: '__MSG_cmdToggleMode__',
      },
      'translate-selection': {
        suggested_key: { default: 'Alt+Z' },
        description: '__MSG_cmdTranslateSelection__',
      },
      'toggle-hover': {
        description: '__MSG_cmdToggleHover__',
      },
    },
    ...(browser === 'firefox'
      ? {
          browser_specific_settings: {
            gecko: { id: 'bioparrot@bioparrot.app', strict_min_version: '115.0' },
          },
        }
      : {}),
  }),
});
