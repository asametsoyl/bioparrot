/* global Zotero, Services, Cc, Ci, APP_SHUTDOWN */
// BioParrot Zotero eklentisi önyükleyicisi (Zotero 7+).
var chromeHandle;

function install() {}

async function startup({ id, version, rootURI }) {
  await Zotero.initializationPromise;
  const aomStartup = Cc['@mozilla.org/addons/addon-manager-startup;1'].getService(Ci.amIAddonManagerStartup);
  const manifestURI = Services.io.newURI(rootURI + 'manifest.json');
  chromeHandle = aomStartup.registerChrome(manifestURI, [['content', 'bioparrot', rootURI + 'content/']]);
  const ctx = { rootURI, pluginID: id, version, Zotero, Services, Cc, Ci };
  ctx._globalThis = ctx;
  Services.scriptloader.loadSubScript(rootURI + 'content/scripts/bioparrot.js', ctx);
  await Zotero.BioParrot.hooks.onStartup();
}

async function onMainWindowLoad({ window }) {
  Zotero.BioParrot?.hooks.onMainWindowLoad(window);
}

async function onMainWindowUnload({ window }) {
  Zotero.BioParrot?.hooks.onMainWindowUnload(window);
}

function shutdown(data, reason) {
  if (reason === APP_SHUTDOWN) return;
  Zotero.BioParrot?.hooks.onShutdown();
  delete Zotero.BioParrot;
  chromeHandle?.destruct();
  chromeHandle = null;
}

function uninstall() {}
