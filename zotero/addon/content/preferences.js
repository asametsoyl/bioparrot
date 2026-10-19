/* global Zotero */
// BioParrot tercih paneli: arayüz, eklentinin ana betiğindeki initPrefsPane ile oluşturulur.
(function init(tries) {
  var root = document.getElementById('bioparrot-prefs-root');
  if (root && Zotero.BioParrot) {
    Zotero.BioParrot.initPrefsPane(document);
    return;
  }
  if (tries < 50) setTimeout(function () { init(tries + 1); }, 100);
})(0);
