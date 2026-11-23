"use strict";

// Two tabs: the speed-up loop and the dbSNP lane. The hash (#dbsnp) selects a tab, so links can point at one.
(function () {
  const tabs = [...document.querySelectorAll('[role="tab"]')];
  function show(id, focus) {
    for (const tab of tabs) {
      const on = tab.getAttribute("aria-controls") === id;
      tab.setAttribute("aria-selected", String(on));
      tab.tabIndex = on ? 0 : -1;
      document.getElementById(tab.getAttribute("aria-controls")).hidden = !on;
      if (on && focus) tab.focus();
    }
    if (location.hash.slice(1) !== id) history.replaceState(null, "", id === "loop" ? location.pathname : "#" + id);
    window.dispatchEvent(new Event("resize"));
  }
  tabs.forEach((tab, i) => {
    tab.addEventListener("click", () => show(tab.getAttribute("aria-controls")));
    tab.addEventListener("keydown", event => {
      const step = {ArrowRight: 1, ArrowLeft: -1}[event.key];
      if (!step) return;
      event.preventDefault();
      show(tabs[(i + step + tabs.length) % tabs.length].getAttribute("aria-controls"), true);
    });
  });
  show(location.hash === "#dbsnp" ? "dbsnp" : "loop");
})();
