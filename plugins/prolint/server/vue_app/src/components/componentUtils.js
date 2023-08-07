//https://www.appsloveworld.com/vuejs/100/76/scroll-to-a-programmatically-selected-row-in-a-v-data-table
export function scrollParentToChild(parent, child) {
  // Where is the parent on page
  var parentRect = parent.getBoundingClientRect();
  // What can you see?
  var parentViewableArea = {
    height: parent.clientHeight,
    width: parent.clientWidth,
  };

  // Where is the child
  var childRect = child.getBoundingClientRect();
  // marco: "dif" as adjustment of scrolling sensitivity and length
  let dif = childRect.top - childRect.bottom;
  // Is the child viewable?
  var isViewable =
    childRect.top >= parentRect.top - dif * 2 &&
    childRect.bottom <= parentRect.top + parentViewableArea.height;
  // if you can't see the child try to scroll parent
  if (!isViewable) {
    // scroll by offset relative to parent
    let scrolling =
      childRect.top +
      parent.scrollTop -
      parentRect.top -
      childRect.height +
      dif;
    //console.log("scrolling...", scrolling, parent, child);
    parent.scrollTop = scrolling;
  }
}
