// node_modules/d3-array/src/ascending.js
function ascending(a, b) {
  return a == null || b == null ? NaN : a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
}

// node_modules/d3-array/src/descending.js
function descending(a, b) {
  return a == null || b == null ? NaN : b < a ? -1 : b > a ? 1 : b >= a ? 0 : NaN;
}

// node_modules/d3-array/src/bisector.js
function bisector(f) {
  let compare1, compare2, delta;
  if (f.length !== 2) {
    compare1 = ascending;
    compare2 = (d, x2) => ascending(f(d), x2);
    delta = (d, x2) => f(d) - x2;
  } else {
    compare1 = f === ascending || f === descending ? f : zero;
    compare2 = f;
    delta = f;
  }
  function left(a, x2, lo = 0, hi = a.length) {
    if (lo < hi) {
      if (compare1(x2, x2) !== 0)
        return hi;
      do {
        const mid = lo + hi >>> 1;
        if (compare2(a[mid], x2) < 0)
          lo = mid + 1;
        else
          hi = mid;
      } while (lo < hi);
    }
    return lo;
  }
  function right(a, x2, lo = 0, hi = a.length) {
    if (lo < hi) {
      if (compare1(x2, x2) !== 0)
        return hi;
      do {
        const mid = lo + hi >>> 1;
        if (compare2(a[mid], x2) <= 0)
          lo = mid + 1;
        else
          hi = mid;
      } while (lo < hi);
    }
    return lo;
  }
  function center(a, x2, lo = 0, hi = a.length) {
    const i = left(a, x2, lo, hi - 1);
    return i > lo && delta(a[i - 1], x2) > -delta(a[i], x2) ? i - 1 : i;
  }
  return { left, center, right };
}
function zero() {
  return 0;
}

// node_modules/d3-array/src/number.js
function number(x2) {
  return x2 === null ? NaN : +x2;
}

// node_modules/d3-array/src/bisect.js
var ascendingBisect = bisector(ascending);
var bisectRight = ascendingBisect.right;
var bisectLeft = ascendingBisect.left;
var bisectCenter = bisector(number).center;
var bisect_default = bisectRight;

// node_modules/d3-array/src/ticks.js
var e10 = Math.sqrt(50);
var e5 = Math.sqrt(10);
var e2 = Math.sqrt(2);
function tickSpec(start2, stop, count) {
  const step = (stop - start2) / Math.max(0, count), power = Math.floor(Math.log10(step)), error = step / Math.pow(10, power), factor = error >= e10 ? 10 : error >= e5 ? 5 : error >= e2 ? 2 : 1;
  let i1, i2, inc;
  if (power < 0) {
    inc = Math.pow(10, -power) / factor;
    i1 = Math.round(start2 * inc);
    i2 = Math.round(stop * inc);
    if (i1 / inc < start2)
      ++i1;
    if (i2 / inc > stop)
      --i2;
    inc = -inc;
  } else {
    inc = Math.pow(10, power) * factor;
    i1 = Math.round(start2 / inc);
    i2 = Math.round(stop / inc);
    if (i1 * inc < start2)
      ++i1;
    if (i2 * inc > stop)
      --i2;
  }
  if (i2 < i1 && 0.5 <= count && count < 2)
    return tickSpec(start2, stop, count * 2);
  return [i1, i2, inc];
}
function ticks(start2, stop, count) {
  stop = +stop, start2 = +start2, count = +count;
  if (!(count > 0))
    return [];
  if (start2 === stop)
    return [start2];
  const reverse = stop < start2, [i1, i2, inc] = reverse ? tickSpec(stop, start2, count) : tickSpec(start2, stop, count);
  if (!(i2 >= i1))
    return [];
  const n = i2 - i1 + 1, ticks2 = new Array(n);
  if (reverse) {
    if (inc < 0)
      for (let i = 0; i < n; ++i)
        ticks2[i] = (i2 - i) / -inc;
    else
      for (let i = 0; i < n; ++i)
        ticks2[i] = (i2 - i) * inc;
  } else {
    if (inc < 0)
      for (let i = 0; i < n; ++i)
        ticks2[i] = (i1 + i) / -inc;
    else
      for (let i = 0; i < n; ++i)
        ticks2[i] = (i1 + i) * inc;
  }
  return ticks2;
}
function tickIncrement(start2, stop, count) {
  stop = +stop, start2 = +start2, count = +count;
  return tickSpec(start2, stop, count)[2];
}
function tickStep(start2, stop, count) {
  stop = +stop, start2 = +start2, count = +count;
  const reverse = stop < start2, inc = reverse ? tickIncrement(stop, start2, count) : tickIncrement(start2, stop, count);
  return (reverse ? -1 : 1) * (inc < 0 ? 1 / -inc : inc);
}

// node_modules/d3-dispatch/src/dispatch.js
var noop = { value: () => {
} };
function dispatch() {
  for (var i = 0, n = arguments.length, _ = {}, t; i < n; ++i) {
    if (!(t = arguments[i] + "") || t in _ || /[\s.]/.test(t))
      throw new Error("illegal type: " + t);
    _[t] = [];
  }
  return new Dispatch(_);
}
function Dispatch(_) {
  this._ = _;
}
function parseTypenames(typenames, types) {
  return typenames.trim().split(/^|\s+/).map(function(t) {
    var name = "", i = t.indexOf(".");
    if (i >= 0)
      name = t.slice(i + 1), t = t.slice(0, i);
    if (t && !types.hasOwnProperty(t))
      throw new Error("unknown type: " + t);
    return { type: t, name };
  });
}
Dispatch.prototype = dispatch.prototype = {
  constructor: Dispatch,
  on: function(typename, callback) {
    var _ = this._, T = parseTypenames(typename + "", _), t, i = -1, n = T.length;
    if (arguments.length < 2) {
      while (++i < n)
        if ((t = (typename = T[i]).type) && (t = get(_[t], typename.name)))
          return t;
      return;
    }
    if (callback != null && typeof callback !== "function")
      throw new Error("invalid callback: " + callback);
    while (++i < n) {
      if (t = (typename = T[i]).type)
        _[t] = set(_[t], typename.name, callback);
      else if (callback == null)
        for (t in _)
          _[t] = set(_[t], typename.name, null);
    }
    return this;
  },
  copy: function() {
    var copy2 = {}, _ = this._;
    for (var t in _)
      copy2[t] = _[t].slice();
    return new Dispatch(copy2);
  },
  call: function(type2, that) {
    if ((n = arguments.length - 2) > 0)
      for (var args = new Array(n), i = 0, n, t; i < n; ++i)
        args[i] = arguments[i + 2];
    if (!this._.hasOwnProperty(type2))
      throw new Error("unknown type: " + type2);
    for (t = this._[type2], i = 0, n = t.length; i < n; ++i)
      t[i].value.apply(that, args);
  },
  apply: function(type2, that, args) {
    if (!this._.hasOwnProperty(type2))
      throw new Error("unknown type: " + type2);
    for (var t = this._[type2], i = 0, n = t.length; i < n; ++i)
      t[i].value.apply(that, args);
  }
};
function get(type2, name) {
  for (var i = 0, n = type2.length, c; i < n; ++i) {
    if ((c = type2[i]).name === name) {
      return c.value;
    }
  }
}
function set(type2, name, callback) {
  for (var i = 0, n = type2.length; i < n; ++i) {
    if (type2[i].name === name) {
      type2[i] = noop, type2 = type2.slice(0, i).concat(type2.slice(i + 1));
      break;
    }
  }
  if (callback != null)
    type2.push({ name, value: callback });
  return type2;
}
var dispatch_default = dispatch;

// node_modules/d3-selection/src/namespaces.js
var xhtml = "http://www.w3.org/1999/xhtml";
var namespaces_default = {
  svg: "http://www.w3.org/2000/svg",
  xhtml,
  xlink: "http://www.w3.org/1999/xlink",
  xml: "http://www.w3.org/XML/1998/namespace",
  xmlns: "http://www.w3.org/2000/xmlns/"
};

// node_modules/d3-selection/src/namespace.js
function namespace_default(name) {
  var prefix = name += "", i = prefix.indexOf(":");
  if (i >= 0 && (prefix = name.slice(0, i)) !== "xmlns")
    name = name.slice(i + 1);
  return namespaces_default.hasOwnProperty(prefix) ? { space: namespaces_default[prefix], local: name } : name;
}

// node_modules/d3-selection/src/creator.js
function creatorInherit(name) {
  return function() {
    var document2 = this.ownerDocument, uri = this.namespaceURI;
    return uri === xhtml && document2.documentElement.namespaceURI === xhtml ? document2.createElement(name) : document2.createElementNS(uri, name);
  };
}
function creatorFixed(fullname) {
  return function() {
    return this.ownerDocument.createElementNS(fullname.space, fullname.local);
  };
}
function creator_default(name) {
  var fullname = namespace_default(name);
  return (fullname.local ? creatorFixed : creatorInherit)(fullname);
}

// node_modules/d3-selection/src/selector.js
function none() {
}
function selector_default(selector) {
  return selector == null ? none : function() {
    return this.querySelector(selector);
  };
}

// node_modules/d3-selection/src/selection/select.js
function select_default(select) {
  if (typeof select !== "function")
    select = selector_default(select);
  for (var groups = this._groups, m = groups.length, subgroups = new Array(m), j = 0; j < m; ++j) {
    for (var group = groups[j], n = group.length, subgroup = subgroups[j] = new Array(n), node, subnode, i = 0; i < n; ++i) {
      if ((node = group[i]) && (subnode = select.call(node, node.__data__, i, group))) {
        if ("__data__" in node)
          subnode.__data__ = node.__data__;
        subgroup[i] = subnode;
      }
    }
  }
  return new Selection(subgroups, this._parents);
}

// node_modules/d3-selection/src/array.js
function array(x2) {
  return x2 == null ? [] : Array.isArray(x2) ? x2 : Array.from(x2);
}

// node_modules/d3-selection/src/selectorAll.js
function empty() {
  return [];
}
function selectorAll_default(selector) {
  return selector == null ? empty : function() {
    return this.querySelectorAll(selector);
  };
}

// node_modules/d3-selection/src/selection/selectAll.js
function arrayAll(select) {
  return function() {
    return array(select.apply(this, arguments));
  };
}
function selectAll_default(select) {
  if (typeof select === "function")
    select = arrayAll(select);
  else
    select = selectorAll_default(select);
  for (var groups = this._groups, m = groups.length, subgroups = [], parents = [], j = 0; j < m; ++j) {
    for (var group = groups[j], n = group.length, node, i = 0; i < n; ++i) {
      if (node = group[i]) {
        subgroups.push(select.call(node, node.__data__, i, group));
        parents.push(node);
      }
    }
  }
  return new Selection(subgroups, parents);
}

// node_modules/d3-selection/src/matcher.js
function matcher_default(selector) {
  return function() {
    return this.matches(selector);
  };
}
function childMatcher(selector) {
  return function(node) {
    return node.matches(selector);
  };
}

// node_modules/d3-selection/src/selection/selectChild.js
var find = Array.prototype.find;
function childFind(match) {
  return function() {
    return find.call(this.children, match);
  };
}
function childFirst() {
  return this.firstElementChild;
}
function selectChild_default(match) {
  return this.select(match == null ? childFirst : childFind(typeof match === "function" ? match : childMatcher(match)));
}

// node_modules/d3-selection/src/selection/selectChildren.js
var filter = Array.prototype.filter;
function children() {
  return Array.from(this.children);
}
function childrenFilter(match) {
  return function() {
    return filter.call(this.children, match);
  };
}
function selectChildren_default(match) {
  return this.selectAll(match == null ? children : childrenFilter(typeof match === "function" ? match : childMatcher(match)));
}

// node_modules/d3-selection/src/selection/filter.js
function filter_default(match) {
  if (typeof match !== "function")
    match = matcher_default(match);
  for (var groups = this._groups, m = groups.length, subgroups = new Array(m), j = 0; j < m; ++j) {
    for (var group = groups[j], n = group.length, subgroup = subgroups[j] = [], node, i = 0; i < n; ++i) {
      if ((node = group[i]) && match.call(node, node.__data__, i, group)) {
        subgroup.push(node);
      }
    }
  }
  return new Selection(subgroups, this._parents);
}

// node_modules/d3-selection/src/selection/sparse.js
function sparse_default(update) {
  return new Array(update.length);
}

// node_modules/d3-selection/src/selection/enter.js
function enter_default() {
  return new Selection(this._enter || this._groups.map(sparse_default), this._parents);
}
function EnterNode(parent, datum2) {
  this.ownerDocument = parent.ownerDocument;
  this.namespaceURI = parent.namespaceURI;
  this._next = null;
  this._parent = parent;
  this.__data__ = datum2;
}
EnterNode.prototype = {
  constructor: EnterNode,
  appendChild: function(child) {
    return this._parent.insertBefore(child, this._next);
  },
  insertBefore: function(child, next) {
    return this._parent.insertBefore(child, next);
  },
  querySelector: function(selector) {
    return this._parent.querySelector(selector);
  },
  querySelectorAll: function(selector) {
    return this._parent.querySelectorAll(selector);
  }
};

// node_modules/d3-selection/src/constant.js
function constant_default(x2) {
  return function() {
    return x2;
  };
}

// node_modules/d3-selection/src/selection/data.js
function bindIndex(parent, group, enter, update, exit, data) {
  var i = 0, node, groupLength = group.length, dataLength = data.length;
  for (; i < dataLength; ++i) {
    if (node = group[i]) {
      node.__data__ = data[i];
      update[i] = node;
    } else {
      enter[i] = new EnterNode(parent, data[i]);
    }
  }
  for (; i < groupLength; ++i) {
    if (node = group[i]) {
      exit[i] = node;
    }
  }
}
function bindKey(parent, group, enter, update, exit, data, key) {
  var i, node, nodeByKeyValue = /* @__PURE__ */ new Map(), groupLength = group.length, dataLength = data.length, keyValues = new Array(groupLength), keyValue;
  for (i = 0; i < groupLength; ++i) {
    if (node = group[i]) {
      keyValues[i] = keyValue = key.call(node, node.__data__, i, group) + "";
      if (nodeByKeyValue.has(keyValue)) {
        exit[i] = node;
      } else {
        nodeByKeyValue.set(keyValue, node);
      }
    }
  }
  for (i = 0; i < dataLength; ++i) {
    keyValue = key.call(parent, data[i], i, data) + "";
    if (node = nodeByKeyValue.get(keyValue)) {
      update[i] = node;
      node.__data__ = data[i];
      nodeByKeyValue.delete(keyValue);
    } else {
      enter[i] = new EnterNode(parent, data[i]);
    }
  }
  for (i = 0; i < groupLength; ++i) {
    if ((node = group[i]) && nodeByKeyValue.get(keyValues[i]) === node) {
      exit[i] = node;
    }
  }
}
function datum(node) {
  return node.__data__;
}
function data_default(value, key) {
  if (!arguments.length)
    return Array.from(this, datum);
  var bind = key ? bindKey : bindIndex, parents = this._parents, groups = this._groups;
  if (typeof value !== "function")
    value = constant_default(value);
  for (var m = groups.length, update = new Array(m), enter = new Array(m), exit = new Array(m), j = 0; j < m; ++j) {
    var parent = parents[j], group = groups[j], groupLength = group.length, data = arraylike(value.call(parent, parent && parent.__data__, j, parents)), dataLength = data.length, enterGroup = enter[j] = new Array(dataLength), updateGroup = update[j] = new Array(dataLength), exitGroup = exit[j] = new Array(groupLength);
    bind(parent, group, enterGroup, updateGroup, exitGroup, data, key);
    for (var i0 = 0, i1 = 0, previous, next; i0 < dataLength; ++i0) {
      if (previous = enterGroup[i0]) {
        if (i0 >= i1)
          i1 = i0 + 1;
        while (!(next = updateGroup[i1]) && ++i1 < dataLength)
          ;
        previous._next = next || null;
      }
    }
  }
  update = new Selection(update, parents);
  update._enter = enter;
  update._exit = exit;
  return update;
}
function arraylike(data) {
  return typeof data === "object" && "length" in data ? data : Array.from(data);
}

// node_modules/d3-selection/src/selection/exit.js
function exit_default() {
  return new Selection(this._exit || this._groups.map(sparse_default), this._parents);
}

// node_modules/d3-selection/src/selection/join.js
function join_default(onenter, onupdate, onexit) {
  var enter = this.enter(), update = this, exit = this.exit();
  if (typeof onenter === "function") {
    enter = onenter(enter);
    if (enter)
      enter = enter.selection();
  } else {
    enter = enter.append(onenter + "");
  }
  if (onupdate != null) {
    update = onupdate(update);
    if (update)
      update = update.selection();
  }
  if (onexit == null)
    exit.remove();
  else
    onexit(exit);
  return enter && update ? enter.merge(update).order() : update;
}

// node_modules/d3-selection/src/selection/merge.js
function merge_default(context) {
  var selection2 = context.selection ? context.selection() : context;
  for (var groups0 = this._groups, groups1 = selection2._groups, m0 = groups0.length, m1 = groups1.length, m = Math.min(m0, m1), merges = new Array(m0), j = 0; j < m; ++j) {
    for (var group0 = groups0[j], group1 = groups1[j], n = group0.length, merge = merges[j] = new Array(n), node, i = 0; i < n; ++i) {
      if (node = group0[i] || group1[i]) {
        merge[i] = node;
      }
    }
  }
  for (; j < m0; ++j) {
    merges[j] = groups0[j];
  }
  return new Selection(merges, this._parents);
}

// node_modules/d3-selection/src/selection/order.js
function order_default() {
  for (var groups = this._groups, j = -1, m = groups.length; ++j < m; ) {
    for (var group = groups[j], i = group.length - 1, next = group[i], node; --i >= 0; ) {
      if (node = group[i]) {
        if (next && node.compareDocumentPosition(next) ^ 4)
          next.parentNode.insertBefore(node, next);
        next = node;
      }
    }
  }
  return this;
}

// node_modules/d3-selection/src/selection/sort.js
function sort_default(compare) {
  if (!compare)
    compare = ascending2;
  function compareNode(a, b) {
    return a && b ? compare(a.__data__, b.__data__) : !a - !b;
  }
  for (var groups = this._groups, m = groups.length, sortgroups = new Array(m), j = 0; j < m; ++j) {
    for (var group = groups[j], n = group.length, sortgroup = sortgroups[j] = new Array(n), node, i = 0; i < n; ++i) {
      if (node = group[i]) {
        sortgroup[i] = node;
      }
    }
    sortgroup.sort(compareNode);
  }
  return new Selection(sortgroups, this._parents).order();
}
function ascending2(a, b) {
  return a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
}

// node_modules/d3-selection/src/selection/call.js
function call_default() {
  var callback = arguments[0];
  arguments[0] = this;
  callback.apply(null, arguments);
  return this;
}

// node_modules/d3-selection/src/selection/nodes.js
function nodes_default() {
  return Array.from(this);
}

// node_modules/d3-selection/src/selection/node.js
function node_default() {
  for (var groups = this._groups, j = 0, m = groups.length; j < m; ++j) {
    for (var group = groups[j], i = 0, n = group.length; i < n; ++i) {
      var node = group[i];
      if (node)
        return node;
    }
  }
  return null;
}

// node_modules/d3-selection/src/selection/size.js
function size_default() {
  let size = 0;
  for (const node of this)
    ++size;
  return size;
}

// node_modules/d3-selection/src/selection/empty.js
function empty_default() {
  return !this.node();
}

// node_modules/d3-selection/src/selection/each.js
function each_default(callback) {
  for (var groups = this._groups, j = 0, m = groups.length; j < m; ++j) {
    for (var group = groups[j], i = 0, n = group.length, node; i < n; ++i) {
      if (node = group[i])
        callback.call(node, node.__data__, i, group);
    }
  }
  return this;
}

// node_modules/d3-selection/src/selection/attr.js
function attrRemove(name) {
  return function() {
    this.removeAttribute(name);
  };
}
function attrRemoveNS(fullname) {
  return function() {
    this.removeAttributeNS(fullname.space, fullname.local);
  };
}
function attrConstant(name, value) {
  return function() {
    this.setAttribute(name, value);
  };
}
function attrConstantNS(fullname, value) {
  return function() {
    this.setAttributeNS(fullname.space, fullname.local, value);
  };
}
function attrFunction(name, value) {
  return function() {
    var v = value.apply(this, arguments);
    if (v == null)
      this.removeAttribute(name);
    else
      this.setAttribute(name, v);
  };
}
function attrFunctionNS(fullname, value) {
  return function() {
    var v = value.apply(this, arguments);
    if (v == null)
      this.removeAttributeNS(fullname.space, fullname.local);
    else
      this.setAttributeNS(fullname.space, fullname.local, v);
  };
}
function attr_default(name, value) {
  var fullname = namespace_default(name);
  if (arguments.length < 2) {
    var node = this.node();
    return fullname.local ? node.getAttributeNS(fullname.space, fullname.local) : node.getAttribute(fullname);
  }
  return this.each((value == null ? fullname.local ? attrRemoveNS : attrRemove : typeof value === "function" ? fullname.local ? attrFunctionNS : attrFunction : fullname.local ? attrConstantNS : attrConstant)(fullname, value));
}

// node_modules/d3-selection/src/window.js
function window_default(node) {
  return node.ownerDocument && node.ownerDocument.defaultView || node.document && node || node.defaultView;
}

// node_modules/d3-selection/src/selection/style.js
function styleRemove(name) {
  return function() {
    this.style.removeProperty(name);
  };
}
function styleConstant(name, value, priority) {
  return function() {
    this.style.setProperty(name, value, priority);
  };
}
function styleFunction(name, value, priority) {
  return function() {
    var v = value.apply(this, arguments);
    if (v == null)
      this.style.removeProperty(name);
    else
      this.style.setProperty(name, v, priority);
  };
}
function style_default(name, value, priority) {
  return arguments.length > 1 ? this.each((value == null ? styleRemove : typeof value === "function" ? styleFunction : styleConstant)(name, value, priority == null ? "" : priority)) : styleValue(this.node(), name);
}
function styleValue(node, name) {
  return node.style.getPropertyValue(name) || window_default(node).getComputedStyle(node, null).getPropertyValue(name);
}

// node_modules/d3-selection/src/selection/property.js
function propertyRemove(name) {
  return function() {
    delete this[name];
  };
}
function propertyConstant(name, value) {
  return function() {
    this[name] = value;
  };
}
function propertyFunction(name, value) {
  return function() {
    var v = value.apply(this, arguments);
    if (v == null)
      delete this[name];
    else
      this[name] = v;
  };
}
function property_default(name, value) {
  return arguments.length > 1 ? this.each((value == null ? propertyRemove : typeof value === "function" ? propertyFunction : propertyConstant)(name, value)) : this.node()[name];
}

// node_modules/d3-selection/src/selection/classed.js
function classArray(string) {
  return string.trim().split(/^|\s+/);
}
function classList(node) {
  return node.classList || new ClassList(node);
}
function ClassList(node) {
  this._node = node;
  this._names = classArray(node.getAttribute("class") || "");
}
ClassList.prototype = {
  add: function(name) {
    var i = this._names.indexOf(name);
    if (i < 0) {
      this._names.push(name);
      this._node.setAttribute("class", this._names.join(" "));
    }
  },
  remove: function(name) {
    var i = this._names.indexOf(name);
    if (i >= 0) {
      this._names.splice(i, 1);
      this._node.setAttribute("class", this._names.join(" "));
    }
  },
  contains: function(name) {
    return this._names.indexOf(name) >= 0;
  }
};
function classedAdd(node, names) {
  var list = classList(node), i = -1, n = names.length;
  while (++i < n)
    list.add(names[i]);
}
function classedRemove(node, names) {
  var list = classList(node), i = -1, n = names.length;
  while (++i < n)
    list.remove(names[i]);
}
function classedTrue(names) {
  return function() {
    classedAdd(this, names);
  };
}
function classedFalse(names) {
  return function() {
    classedRemove(this, names);
  };
}
function classedFunction(names, value) {
  return function() {
    (value.apply(this, arguments) ? classedAdd : classedRemove)(this, names);
  };
}
function classed_default(name, value) {
  var names = classArray(name + "");
  if (arguments.length < 2) {
    var list = classList(this.node()), i = -1, n = names.length;
    while (++i < n)
      if (!list.contains(names[i]))
        return false;
    return true;
  }
  return this.each((typeof value === "function" ? classedFunction : value ? classedTrue : classedFalse)(names, value));
}

// node_modules/d3-selection/src/selection/text.js
function textRemove() {
  this.textContent = "";
}
function textConstant(value) {
  return function() {
    this.textContent = value;
  };
}
function textFunction(value) {
  return function() {
    var v = value.apply(this, arguments);
    this.textContent = v == null ? "" : v;
  };
}
function text_default(value) {
  return arguments.length ? this.each(value == null ? textRemove : (typeof value === "function" ? textFunction : textConstant)(value)) : this.node().textContent;
}

// node_modules/d3-selection/src/selection/html.js
function htmlRemove() {
  this.innerHTML = "";
}
function htmlConstant(value) {
  return function() {
    this.innerHTML = value;
  };
}
function htmlFunction(value) {
  return function() {
    var v = value.apply(this, arguments);
    this.innerHTML = v == null ? "" : v;
  };
}
function html_default(value) {
  return arguments.length ? this.each(value == null ? htmlRemove : (typeof value === "function" ? htmlFunction : htmlConstant)(value)) : this.node().innerHTML;
}

// node_modules/d3-selection/src/selection/raise.js
function raise() {
  if (this.nextSibling)
    this.parentNode.appendChild(this);
}
function raise_default() {
  return this.each(raise);
}

// node_modules/d3-selection/src/selection/lower.js
function lower() {
  if (this.previousSibling)
    this.parentNode.insertBefore(this, this.parentNode.firstChild);
}
function lower_default() {
  return this.each(lower);
}

// node_modules/d3-selection/src/selection/append.js
function append_default(name) {
  var create2 = typeof name === "function" ? name : creator_default(name);
  return this.select(function() {
    return this.appendChild(create2.apply(this, arguments));
  });
}

// node_modules/d3-selection/src/selection/insert.js
function constantNull() {
  return null;
}
function insert_default(name, before) {
  var create2 = typeof name === "function" ? name : creator_default(name), select = before == null ? constantNull : typeof before === "function" ? before : selector_default(before);
  return this.select(function() {
    return this.insertBefore(create2.apply(this, arguments), select.apply(this, arguments) || null);
  });
}

// node_modules/d3-selection/src/selection/remove.js
function remove() {
  var parent = this.parentNode;
  if (parent)
    parent.removeChild(this);
}
function remove_default() {
  return this.each(remove);
}

// node_modules/d3-selection/src/selection/clone.js
function selection_cloneShallow() {
  var clone = this.cloneNode(false), parent = this.parentNode;
  return parent ? parent.insertBefore(clone, this.nextSibling) : clone;
}
function selection_cloneDeep() {
  var clone = this.cloneNode(true), parent = this.parentNode;
  return parent ? parent.insertBefore(clone, this.nextSibling) : clone;
}
function clone_default(deep) {
  return this.select(deep ? selection_cloneDeep : selection_cloneShallow);
}

// node_modules/d3-selection/src/selection/datum.js
function datum_default(value) {
  return arguments.length ? this.property("__data__", value) : this.node().__data__;
}

// node_modules/d3-selection/src/selection/on.js
function contextListener(listener) {
  return function(event) {
    listener.call(this, event, this.__data__);
  };
}
function parseTypenames2(typenames) {
  return typenames.trim().split(/^|\s+/).map(function(t) {
    var name = "", i = t.indexOf(".");
    if (i >= 0)
      name = t.slice(i + 1), t = t.slice(0, i);
    return { type: t, name };
  });
}
function onRemove(typename) {
  return function() {
    var on = this.__on;
    if (!on)
      return;
    for (var j = 0, i = -1, m = on.length, o; j < m; ++j) {
      if (o = on[j], (!typename.type || o.type === typename.type) && o.name === typename.name) {
        this.removeEventListener(o.type, o.listener, o.options);
      } else {
        on[++i] = o;
      }
    }
    if (++i)
      on.length = i;
    else
      delete this.__on;
  };
}
function onAdd(typename, value, options) {
  return function() {
    var on = this.__on, o, listener = contextListener(value);
    if (on)
      for (var j = 0, m = on.length; j < m; ++j) {
        if ((o = on[j]).type === typename.type && o.name === typename.name) {
          this.removeEventListener(o.type, o.listener, o.options);
          this.addEventListener(o.type, o.listener = listener, o.options = options);
          o.value = value;
          return;
        }
      }
    this.addEventListener(typename.type, listener, options);
    o = { type: typename.type, name: typename.name, value, listener, options };
    if (!on)
      this.__on = [o];
    else
      on.push(o);
  };
}
function on_default(typename, value, options) {
  var typenames = parseTypenames2(typename + ""), i, n = typenames.length, t;
  if (arguments.length < 2) {
    var on = this.node().__on;
    if (on)
      for (var j = 0, m = on.length, o; j < m; ++j) {
        for (i = 0, o = on[j]; i < n; ++i) {
          if ((t = typenames[i]).type === o.type && t.name === o.name) {
            return o.value;
          }
        }
      }
    return;
  }
  on = value ? onAdd : onRemove;
  for (i = 0; i < n; ++i)
    this.each(on(typenames[i], value, options));
  return this;
}

// node_modules/d3-selection/src/selection/dispatch.js
function dispatchEvent(node, type2, params) {
  var window2 = window_default(node), event = window2.CustomEvent;
  if (typeof event === "function") {
    event = new event(type2, params);
  } else {
    event = window2.document.createEvent("Event");
    if (params)
      event.initEvent(type2, params.bubbles, params.cancelable), event.detail = params.detail;
    else
      event.initEvent(type2, false, false);
  }
  node.dispatchEvent(event);
}
function dispatchConstant(type2, params) {
  return function() {
    return dispatchEvent(this, type2, params);
  };
}
function dispatchFunction(type2, params) {
  return function() {
    return dispatchEvent(this, type2, params.apply(this, arguments));
  };
}
function dispatch_default2(type2, params) {
  return this.each((typeof params === "function" ? dispatchFunction : dispatchConstant)(type2, params));
}

// node_modules/d3-selection/src/selection/iterator.js
function* iterator_default() {
  for (var groups = this._groups, j = 0, m = groups.length; j < m; ++j) {
    for (var group = groups[j], i = 0, n = group.length, node; i < n; ++i) {
      if (node = group[i])
        yield node;
    }
  }
}

// node_modules/d3-selection/src/selection/index.js
var root = [null];
function Selection(groups, parents) {
  this._groups = groups;
  this._parents = parents;
}
function selection() {
  return new Selection([[document.documentElement]], root);
}
function selection_selection() {
  return this;
}
Selection.prototype = selection.prototype = {
  constructor: Selection,
  select: select_default,
  selectAll: selectAll_default,
  selectChild: selectChild_default,
  selectChildren: selectChildren_default,
  filter: filter_default,
  data: data_default,
  enter: enter_default,
  exit: exit_default,
  join: join_default,
  merge: merge_default,
  selection: selection_selection,
  order: order_default,
  sort: sort_default,
  call: call_default,
  nodes: nodes_default,
  node: node_default,
  size: size_default,
  empty: empty_default,
  each: each_default,
  attr: attr_default,
  style: style_default,
  property: property_default,
  classed: classed_default,
  text: text_default,
  html: html_default,
  raise: raise_default,
  lower: lower_default,
  append: append_default,
  insert: insert_default,
  remove: remove_default,
  clone: clone_default,
  datum: datum_default,
  on: on_default,
  dispatch: dispatch_default2,
  [Symbol.iterator]: iterator_default
};
var selection_default = selection;

// node_modules/d3-selection/src/select.js
function select_default2(selector) {
  return typeof selector === "string" ? new Selection([[document.querySelector(selector)]], [document.documentElement]) : new Selection([[selector]], root);
}

// node_modules/d3-selection/src/sourceEvent.js
function sourceEvent_default(event) {
  let sourceEvent;
  while (sourceEvent = event.sourceEvent)
    event = sourceEvent;
  return event;
}

// node_modules/d3-selection/src/pointer.js
function pointer_default(event, node) {
  event = sourceEvent_default(event);
  if (node === void 0)
    node = event.currentTarget;
  if (node) {
    var svg = node.ownerSVGElement || node;
    if (svg.createSVGPoint) {
      var point3 = svg.createSVGPoint();
      point3.x = event.clientX, point3.y = event.clientY;
      point3 = point3.matrixTransform(node.getScreenCTM().inverse());
      return [point3.x, point3.y];
    }
    if (node.getBoundingClientRect) {
      var rect = node.getBoundingClientRect();
      return [event.clientX - rect.left - node.clientLeft, event.clientY - rect.top - node.clientTop];
    }
  }
  return [event.pageX, event.pageY];
}

// node_modules/d3-drag/src/noevent.js
var nonpassivecapture = { capture: true, passive: false };
function noevent_default(event) {
  event.preventDefault();
  event.stopImmediatePropagation();
}

// node_modules/d3-drag/src/nodrag.js
function nodrag_default(view) {
  var root2 = view.document.documentElement, selection2 = select_default2(view).on("dragstart.drag", noevent_default, nonpassivecapture);
  if ("onselectstart" in root2) {
    selection2.on("selectstart.drag", noevent_default, nonpassivecapture);
  } else {
    root2.__noselect = root2.style.MozUserSelect;
    root2.style.MozUserSelect = "none";
  }
}
function yesdrag(view, noclick) {
  var root2 = view.document.documentElement, selection2 = select_default2(view).on("dragstart.drag", null);
  if (noclick) {
    selection2.on("click.drag", noevent_default, nonpassivecapture);
    setTimeout(function() {
      selection2.on("click.drag", null);
    }, 0);
  }
  if ("onselectstart" in root2) {
    selection2.on("selectstart.drag", null);
  } else {
    root2.style.MozUserSelect = root2.__noselect;
    delete root2.__noselect;
  }
}

// node_modules/d3-color/src/define.js
function define_default(constructor, factory, prototype) {
  constructor.prototype = factory.prototype = prototype;
  prototype.constructor = constructor;
}
function extend(parent, definition) {
  var prototype = Object.create(parent.prototype);
  for (var key in definition)
    prototype[key] = definition[key];
  return prototype;
}

// node_modules/d3-color/src/color.js
function Color() {
}
var darker = 0.7;
var brighter = 1 / darker;
var reI = "\\s*([+-]?\\d+)\\s*";
var reN = "\\s*([+-]?(?:\\d*\\.)?\\d+(?:[eE][+-]?\\d+)?)\\s*";
var reP = "\\s*([+-]?(?:\\d*\\.)?\\d+(?:[eE][+-]?\\d+)?)%\\s*";
var reHex = /^#([0-9a-f]{3,8})$/;
var reRgbInteger = new RegExp(`^rgb\\(${reI},${reI},${reI}\\)$`);
var reRgbPercent = new RegExp(`^rgb\\(${reP},${reP},${reP}\\)$`);
var reRgbaInteger = new RegExp(`^rgba\\(${reI},${reI},${reI},${reN}\\)$`);
var reRgbaPercent = new RegExp(`^rgba\\(${reP},${reP},${reP},${reN}\\)$`);
var reHslPercent = new RegExp(`^hsl\\(${reN},${reP},${reP}\\)$`);
var reHslaPercent = new RegExp(`^hsla\\(${reN},${reP},${reP},${reN}\\)$`);
var named = {
  aliceblue: 15792383,
  antiquewhite: 16444375,
  aqua: 65535,
  aquamarine: 8388564,
  azure: 15794175,
  beige: 16119260,
  bisque: 16770244,
  black: 0,
  blanchedalmond: 16772045,
  blue: 255,
  blueviolet: 9055202,
  brown: 10824234,
  burlywood: 14596231,
  cadetblue: 6266528,
  chartreuse: 8388352,
  chocolate: 13789470,
  coral: 16744272,
  cornflowerblue: 6591981,
  cornsilk: 16775388,
  crimson: 14423100,
  cyan: 65535,
  darkblue: 139,
  darkcyan: 35723,
  darkgoldenrod: 12092939,
  darkgray: 11119017,
  darkgreen: 25600,
  darkgrey: 11119017,
  darkkhaki: 12433259,
  darkmagenta: 9109643,
  darkolivegreen: 5597999,
  darkorange: 16747520,
  darkorchid: 10040012,
  darkred: 9109504,
  darksalmon: 15308410,
  darkseagreen: 9419919,
  darkslateblue: 4734347,
  darkslategray: 3100495,
  darkslategrey: 3100495,
  darkturquoise: 52945,
  darkviolet: 9699539,
  deeppink: 16716947,
  deepskyblue: 49151,
  dimgray: 6908265,
  dimgrey: 6908265,
  dodgerblue: 2003199,
  firebrick: 11674146,
  floralwhite: 16775920,
  forestgreen: 2263842,
  fuchsia: 16711935,
  gainsboro: 14474460,
  ghostwhite: 16316671,
  gold: 16766720,
  goldenrod: 14329120,
  gray: 8421504,
  green: 32768,
  greenyellow: 11403055,
  grey: 8421504,
  honeydew: 15794160,
  hotpink: 16738740,
  indianred: 13458524,
  indigo: 4915330,
  ivory: 16777200,
  khaki: 15787660,
  lavender: 15132410,
  lavenderblush: 16773365,
  lawngreen: 8190976,
  lemonchiffon: 16775885,
  lightblue: 11393254,
  lightcoral: 15761536,
  lightcyan: 14745599,
  lightgoldenrodyellow: 16448210,
  lightgray: 13882323,
  lightgreen: 9498256,
  lightgrey: 13882323,
  lightpink: 16758465,
  lightsalmon: 16752762,
  lightseagreen: 2142890,
  lightskyblue: 8900346,
  lightslategray: 7833753,
  lightslategrey: 7833753,
  lightsteelblue: 11584734,
  lightyellow: 16777184,
  lime: 65280,
  limegreen: 3329330,
  linen: 16445670,
  magenta: 16711935,
  maroon: 8388608,
  mediumaquamarine: 6737322,
  mediumblue: 205,
  mediumorchid: 12211667,
  mediumpurple: 9662683,
  mediumseagreen: 3978097,
  mediumslateblue: 8087790,
  mediumspringgreen: 64154,
  mediumturquoise: 4772300,
  mediumvioletred: 13047173,
  midnightblue: 1644912,
  mintcream: 16121850,
  mistyrose: 16770273,
  moccasin: 16770229,
  navajowhite: 16768685,
  navy: 128,
  oldlace: 16643558,
  olive: 8421376,
  olivedrab: 7048739,
  orange: 16753920,
  orangered: 16729344,
  orchid: 14315734,
  palegoldenrod: 15657130,
  palegreen: 10025880,
  paleturquoise: 11529966,
  palevioletred: 14381203,
  papayawhip: 16773077,
  peachpuff: 16767673,
  peru: 13468991,
  pink: 16761035,
  plum: 14524637,
  powderblue: 11591910,
  purple: 8388736,
  rebeccapurple: 6697881,
  red: 16711680,
  rosybrown: 12357519,
  royalblue: 4286945,
  saddlebrown: 9127187,
  salmon: 16416882,
  sandybrown: 16032864,
  seagreen: 3050327,
  seashell: 16774638,
  sienna: 10506797,
  silver: 12632256,
  skyblue: 8900331,
  slateblue: 6970061,
  slategray: 7372944,
  slategrey: 7372944,
  snow: 16775930,
  springgreen: 65407,
  steelblue: 4620980,
  tan: 13808780,
  teal: 32896,
  thistle: 14204888,
  tomato: 16737095,
  turquoise: 4251856,
  violet: 15631086,
  wheat: 16113331,
  white: 16777215,
  whitesmoke: 16119285,
  yellow: 16776960,
  yellowgreen: 10145074
};
define_default(Color, color, {
  copy(channels) {
    return Object.assign(new this.constructor(), this, channels);
  },
  displayable() {
    return this.rgb().displayable();
  },
  hex: color_formatHex,
  // Deprecated! Use color.formatHex.
  formatHex: color_formatHex,
  formatHex8: color_formatHex8,
  formatHsl: color_formatHsl,
  formatRgb: color_formatRgb,
  toString: color_formatRgb
});
function color_formatHex() {
  return this.rgb().formatHex();
}
function color_formatHex8() {
  return this.rgb().formatHex8();
}
function color_formatHsl() {
  return hslConvert(this).formatHsl();
}
function color_formatRgb() {
  return this.rgb().formatRgb();
}
function color(format2) {
  var m, l;
  format2 = (format2 + "").trim().toLowerCase();
  return (m = reHex.exec(format2)) ? (l = m[1].length, m = parseInt(m[1], 16), l === 6 ? rgbn(m) : l === 3 ? new Rgb(m >> 8 & 15 | m >> 4 & 240, m >> 4 & 15 | m & 240, (m & 15) << 4 | m & 15, 1) : l === 8 ? rgba(m >> 24 & 255, m >> 16 & 255, m >> 8 & 255, (m & 255) / 255) : l === 4 ? rgba(m >> 12 & 15 | m >> 8 & 240, m >> 8 & 15 | m >> 4 & 240, m >> 4 & 15 | m & 240, ((m & 15) << 4 | m & 15) / 255) : null) : (m = reRgbInteger.exec(format2)) ? new Rgb(m[1], m[2], m[3], 1) : (m = reRgbPercent.exec(format2)) ? new Rgb(m[1] * 255 / 100, m[2] * 255 / 100, m[3] * 255 / 100, 1) : (m = reRgbaInteger.exec(format2)) ? rgba(m[1], m[2], m[3], m[4]) : (m = reRgbaPercent.exec(format2)) ? rgba(m[1] * 255 / 100, m[2] * 255 / 100, m[3] * 255 / 100, m[4]) : (m = reHslPercent.exec(format2)) ? hsla(m[1], m[2] / 100, m[3] / 100, 1) : (m = reHslaPercent.exec(format2)) ? hsla(m[1], m[2] / 100, m[3] / 100, m[4]) : named.hasOwnProperty(format2) ? rgbn(named[format2]) : format2 === "transparent" ? new Rgb(NaN, NaN, NaN, 0) : null;
}
function rgbn(n) {
  return new Rgb(n >> 16 & 255, n >> 8 & 255, n & 255, 1);
}
function rgba(r, g, b, a) {
  if (a <= 0)
    r = g = b = NaN;
  return new Rgb(r, g, b, a);
}
function rgbConvert(o) {
  if (!(o instanceof Color))
    o = color(o);
  if (!o)
    return new Rgb();
  o = o.rgb();
  return new Rgb(o.r, o.g, o.b, o.opacity);
}
function rgb(r, g, b, opacity) {
  return arguments.length === 1 ? rgbConvert(r) : new Rgb(r, g, b, opacity == null ? 1 : opacity);
}
function Rgb(r, g, b, opacity) {
  this.r = +r;
  this.g = +g;
  this.b = +b;
  this.opacity = +opacity;
}
define_default(Rgb, rgb, extend(Color, {
  brighter(k) {
    k = k == null ? brighter : Math.pow(brighter, k);
    return new Rgb(this.r * k, this.g * k, this.b * k, this.opacity);
  },
  darker(k) {
    k = k == null ? darker : Math.pow(darker, k);
    return new Rgb(this.r * k, this.g * k, this.b * k, this.opacity);
  },
  rgb() {
    return this;
  },
  clamp() {
    return new Rgb(clampi(this.r), clampi(this.g), clampi(this.b), clampa(this.opacity));
  },
  displayable() {
    return -0.5 <= this.r && this.r < 255.5 && (-0.5 <= this.g && this.g < 255.5) && (-0.5 <= this.b && this.b < 255.5) && (0 <= this.opacity && this.opacity <= 1);
  },
  hex: rgb_formatHex,
  // Deprecated! Use color.formatHex.
  formatHex: rgb_formatHex,
  formatHex8: rgb_formatHex8,
  formatRgb: rgb_formatRgb,
  toString: rgb_formatRgb
}));
function rgb_formatHex() {
  return `#${hex(this.r)}${hex(this.g)}${hex(this.b)}`;
}
function rgb_formatHex8() {
  return `#${hex(this.r)}${hex(this.g)}${hex(this.b)}${hex((isNaN(this.opacity) ? 1 : this.opacity) * 255)}`;
}
function rgb_formatRgb() {
  const a = clampa(this.opacity);
  return `${a === 1 ? "rgb(" : "rgba("}${clampi(this.r)}, ${clampi(this.g)}, ${clampi(this.b)}${a === 1 ? ")" : `, ${a})`}`;
}
function clampa(opacity) {
  return isNaN(opacity) ? 1 : Math.max(0, Math.min(1, opacity));
}
function clampi(value) {
  return Math.max(0, Math.min(255, Math.round(value) || 0));
}
function hex(value) {
  value = clampi(value);
  return (value < 16 ? "0" : "") + value.toString(16);
}
function hsla(h, s, l, a) {
  if (a <= 0)
    h = s = l = NaN;
  else if (l <= 0 || l >= 1)
    h = s = NaN;
  else if (s <= 0)
    h = NaN;
  return new Hsl(h, s, l, a);
}
function hslConvert(o) {
  if (o instanceof Hsl)
    return new Hsl(o.h, o.s, o.l, o.opacity);
  if (!(o instanceof Color))
    o = color(o);
  if (!o)
    return new Hsl();
  if (o instanceof Hsl)
    return o;
  o = o.rgb();
  var r = o.r / 255, g = o.g / 255, b = o.b / 255, min2 = Math.min(r, g, b), max2 = Math.max(r, g, b), h = NaN, s = max2 - min2, l = (max2 + min2) / 2;
  if (s) {
    if (r === max2)
      h = (g - b) / s + (g < b) * 6;
    else if (g === max2)
      h = (b - r) / s + 2;
    else
      h = (r - g) / s + 4;
    s /= l < 0.5 ? max2 + min2 : 2 - max2 - min2;
    h *= 60;
  } else {
    s = l > 0 && l < 1 ? 0 : h;
  }
  return new Hsl(h, s, l, o.opacity);
}
function hsl(h, s, l, opacity) {
  return arguments.length === 1 ? hslConvert(h) : new Hsl(h, s, l, opacity == null ? 1 : opacity);
}
function Hsl(h, s, l, opacity) {
  this.h = +h;
  this.s = +s;
  this.l = +l;
  this.opacity = +opacity;
}
define_default(Hsl, hsl, extend(Color, {
  brighter(k) {
    k = k == null ? brighter : Math.pow(brighter, k);
    return new Hsl(this.h, this.s, this.l * k, this.opacity);
  },
  darker(k) {
    k = k == null ? darker : Math.pow(darker, k);
    return new Hsl(this.h, this.s, this.l * k, this.opacity);
  },
  rgb() {
    var h = this.h % 360 + (this.h < 0) * 360, s = isNaN(h) || isNaN(this.s) ? 0 : this.s, l = this.l, m2 = l + (l < 0.5 ? l : 1 - l) * s, m1 = 2 * l - m2;
    return new Rgb(
      hsl2rgb(h >= 240 ? h - 240 : h + 120, m1, m2),
      hsl2rgb(h, m1, m2),
      hsl2rgb(h < 120 ? h + 240 : h - 120, m1, m2),
      this.opacity
    );
  },
  clamp() {
    return new Hsl(clamph(this.h), clampt(this.s), clampt(this.l), clampa(this.opacity));
  },
  displayable() {
    return (0 <= this.s && this.s <= 1 || isNaN(this.s)) && (0 <= this.l && this.l <= 1) && (0 <= this.opacity && this.opacity <= 1);
  },
  formatHsl() {
    const a = clampa(this.opacity);
    return `${a === 1 ? "hsl(" : "hsla("}${clamph(this.h)}, ${clampt(this.s) * 100}%, ${clampt(this.l) * 100}%${a === 1 ? ")" : `, ${a})`}`;
  }
}));
function clamph(value) {
  value = (value || 0) % 360;
  return value < 0 ? value + 360 : value;
}
function clampt(value) {
  return Math.max(0, Math.min(1, value || 0));
}
function hsl2rgb(h, m1, m2) {
  return (h < 60 ? m1 + (m2 - m1) * h / 60 : h < 180 ? m2 : h < 240 ? m1 + (m2 - m1) * (240 - h) / 60 : m1) * 255;
}

// node_modules/d3-interpolate/src/basis.js
function basis(t1, v0, v1, v2, v3) {
  var t2 = t1 * t1, t3 = t2 * t1;
  return ((1 - 3 * t1 + 3 * t2 - t3) * v0 + (4 - 6 * t2 + 3 * t3) * v1 + (1 + 3 * t1 + 3 * t2 - 3 * t3) * v2 + t3 * v3) / 6;
}
function basis_default(values) {
  var n = values.length - 1;
  return function(t) {
    var i = t <= 0 ? t = 0 : t >= 1 ? (t = 1, n - 1) : Math.floor(t * n), v1 = values[i], v2 = values[i + 1], v0 = i > 0 ? values[i - 1] : 2 * v1 - v2, v3 = i < n - 1 ? values[i + 2] : 2 * v2 - v1;
    return basis((t - i / n) * n, v0, v1, v2, v3);
  };
}

// node_modules/d3-interpolate/src/basisClosed.js
function basisClosed_default(values) {
  var n = values.length;
  return function(t) {
    var i = Math.floor(((t %= 1) < 0 ? ++t : t) * n), v0 = values[(i + n - 1) % n], v1 = values[i % n], v2 = values[(i + 1) % n], v3 = values[(i + 2) % n];
    return basis((t - i / n) * n, v0, v1, v2, v3);
  };
}

// node_modules/d3-interpolate/src/constant.js
var constant_default2 = (x2) => () => x2;

// node_modules/d3-interpolate/src/color.js
function linear(a, d) {
  return function(t) {
    return a + t * d;
  };
}
function exponential(a, b, y2) {
  return a = Math.pow(a, y2), b = Math.pow(b, y2) - a, y2 = 1 / y2, function(t) {
    return Math.pow(a + t * b, y2);
  };
}
function gamma(y2) {
  return (y2 = +y2) === 1 ? nogamma : function(a, b) {
    return b - a ? exponential(a, b, y2) : constant_default2(isNaN(a) ? b : a);
  };
}
function nogamma(a, b) {
  var d = b - a;
  return d ? linear(a, d) : constant_default2(isNaN(a) ? b : a);
}

// node_modules/d3-interpolate/src/rgb.js
var rgb_default = function rgbGamma(y2) {
  var color2 = gamma(y2);
  function rgb2(start2, end) {
    var r = color2((start2 = rgb(start2)).r, (end = rgb(end)).r), g = color2(start2.g, end.g), b = color2(start2.b, end.b), opacity = nogamma(start2.opacity, end.opacity);
    return function(t) {
      start2.r = r(t);
      start2.g = g(t);
      start2.b = b(t);
      start2.opacity = opacity(t);
      return start2 + "";
    };
  }
  rgb2.gamma = rgbGamma;
  return rgb2;
}(1);
function rgbSpline(spline) {
  return function(colors) {
    var n = colors.length, r = new Array(n), g = new Array(n), b = new Array(n), i, color2;
    for (i = 0; i < n; ++i) {
      color2 = rgb(colors[i]);
      r[i] = color2.r || 0;
      g[i] = color2.g || 0;
      b[i] = color2.b || 0;
    }
    r = spline(r);
    g = spline(g);
    b = spline(b);
    color2.opacity = 1;
    return function(t) {
      color2.r = r(t);
      color2.g = g(t);
      color2.b = b(t);
      return color2 + "";
    };
  };
}
var rgbBasis = rgbSpline(basis_default);
var rgbBasisClosed = rgbSpline(basisClosed_default);

// node_modules/d3-interpolate/src/numberArray.js
function numberArray_default(a, b) {
  if (!b)
    b = [];
  var n = a ? Math.min(b.length, a.length) : 0, c = b.slice(), i;
  return function(t) {
    for (i = 0; i < n; ++i)
      c[i] = a[i] * (1 - t) + b[i] * t;
    return c;
  };
}
function isNumberArray(x2) {
  return ArrayBuffer.isView(x2) && !(x2 instanceof DataView);
}

// node_modules/d3-interpolate/src/array.js
function genericArray(a, b) {
  var nb = b ? b.length : 0, na = a ? Math.min(nb, a.length) : 0, x2 = new Array(na), c = new Array(nb), i;
  for (i = 0; i < na; ++i)
    x2[i] = value_default(a[i], b[i]);
  for (; i < nb; ++i)
    c[i] = b[i];
  return function(t) {
    for (i = 0; i < na; ++i)
      c[i] = x2[i](t);
    return c;
  };
}

// node_modules/d3-interpolate/src/date.js
function date_default(a, b) {
  var d = /* @__PURE__ */ new Date();
  return a = +a, b = +b, function(t) {
    return d.setTime(a * (1 - t) + b * t), d;
  };
}

// node_modules/d3-interpolate/src/number.js
function number_default(a, b) {
  return a = +a, b = +b, function(t) {
    return a * (1 - t) + b * t;
  };
}

// node_modules/d3-interpolate/src/object.js
function object_default(a, b) {
  var i = {}, c = {}, k;
  if (a === null || typeof a !== "object")
    a = {};
  if (b === null || typeof b !== "object")
    b = {};
  for (k in b) {
    if (k in a) {
      i[k] = value_default(a[k], b[k]);
    } else {
      c[k] = b[k];
    }
  }
  return function(t) {
    for (k in i)
      c[k] = i[k](t);
    return c;
  };
}

// node_modules/d3-interpolate/src/string.js
var reA = /[-+]?(?:\d+\.?\d*|\.?\d+)(?:[eE][-+]?\d+)?/g;
var reB = new RegExp(reA.source, "g");
function zero2(b) {
  return function() {
    return b;
  };
}
function one(b) {
  return function(t) {
    return b(t) + "";
  };
}
function string_default(a, b) {
  var bi = reA.lastIndex = reB.lastIndex = 0, am, bm, bs, i = -1, s = [], q = [];
  a = a + "", b = b + "";
  while ((am = reA.exec(a)) && (bm = reB.exec(b))) {
    if ((bs = bm.index) > bi) {
      bs = b.slice(bi, bs);
      if (s[i])
        s[i] += bs;
      else
        s[++i] = bs;
    }
    if ((am = am[0]) === (bm = bm[0])) {
      if (s[i])
        s[i] += bm;
      else
        s[++i] = bm;
    } else {
      s[++i] = null;
      q.push({ i, x: number_default(am, bm) });
    }
    bi = reB.lastIndex;
  }
  if (bi < b.length) {
    bs = b.slice(bi);
    if (s[i])
      s[i] += bs;
    else
      s[++i] = bs;
  }
  return s.length < 2 ? q[0] ? one(q[0].x) : zero2(b) : (b = q.length, function(t) {
    for (var i2 = 0, o; i2 < b; ++i2)
      s[(o = q[i2]).i] = o.x(t);
    return s.join("");
  });
}

// node_modules/d3-interpolate/src/value.js
function value_default(a, b) {
  var t = typeof b, c;
  return b == null || t === "boolean" ? constant_default2(b) : (t === "number" ? number_default : t === "string" ? (c = color(b)) ? (b = c, rgb_default) : string_default : b instanceof color ? rgb_default : b instanceof Date ? date_default : isNumberArray(b) ? numberArray_default : Array.isArray(b) ? genericArray : typeof b.valueOf !== "function" && typeof b.toString !== "function" || isNaN(b) ? object_default : number_default)(a, b);
}

// node_modules/d3-interpolate/src/round.js
function round_default(a, b) {
  return a = +a, b = +b, function(t) {
    return Math.round(a * (1 - t) + b * t);
  };
}

// node_modules/d3-interpolate/src/transform/decompose.js
var degrees = 180 / Math.PI;
var identity = {
  translateX: 0,
  translateY: 0,
  rotate: 0,
  skewX: 0,
  scaleX: 1,
  scaleY: 1
};
function decompose_default(a, b, c, d, e, f) {
  var scaleX, scaleY, skewX;
  if (scaleX = Math.sqrt(a * a + b * b))
    a /= scaleX, b /= scaleX;
  if (skewX = a * c + b * d)
    c -= a * skewX, d -= b * skewX;
  if (scaleY = Math.sqrt(c * c + d * d))
    c /= scaleY, d /= scaleY, skewX /= scaleY;
  if (a * d < b * c)
    a = -a, b = -b, skewX = -skewX, scaleX = -scaleX;
  return {
    translateX: e,
    translateY: f,
    rotate: Math.atan2(b, a) * degrees,
    skewX: Math.atan(skewX) * degrees,
    scaleX,
    scaleY
  };
}

// node_modules/d3-interpolate/src/transform/parse.js
var svgNode;
function parseCss(value) {
  const m = new (typeof DOMMatrix === "function" ? DOMMatrix : WebKitCSSMatrix)(value + "");
  return m.isIdentity ? identity : decompose_default(m.a, m.b, m.c, m.d, m.e, m.f);
}
function parseSvg(value) {
  if (value == null)
    return identity;
  if (!svgNode)
    svgNode = document.createElementNS("http://www.w3.org/2000/svg", "g");
  svgNode.setAttribute("transform", value);
  if (!(value = svgNode.transform.baseVal.consolidate()))
    return identity;
  value = value.matrix;
  return decompose_default(value.a, value.b, value.c, value.d, value.e, value.f);
}

// node_modules/d3-interpolate/src/transform/index.js
function interpolateTransform(parse, pxComma, pxParen, degParen) {
  function pop(s) {
    return s.length ? s.pop() + " " : "";
  }
  function translate(xa, ya, xb, yb, s, q) {
    if (xa !== xb || ya !== yb) {
      var i = s.push("translate(", null, pxComma, null, pxParen);
      q.push({ i: i - 4, x: number_default(xa, xb) }, { i: i - 2, x: number_default(ya, yb) });
    } else if (xb || yb) {
      s.push("translate(" + xb + pxComma + yb + pxParen);
    }
  }
  function rotate(a, b, s, q) {
    if (a !== b) {
      if (a - b > 180)
        b += 360;
      else if (b - a > 180)
        a += 360;
      q.push({ i: s.push(pop(s) + "rotate(", null, degParen) - 2, x: number_default(a, b) });
    } else if (b) {
      s.push(pop(s) + "rotate(" + b + degParen);
    }
  }
  function skewX(a, b, s, q) {
    if (a !== b) {
      q.push({ i: s.push(pop(s) + "skewX(", null, degParen) - 2, x: number_default(a, b) });
    } else if (b) {
      s.push(pop(s) + "skewX(" + b + degParen);
    }
  }
  function scale(xa, ya, xb, yb, s, q) {
    if (xa !== xb || ya !== yb) {
      var i = s.push(pop(s) + "scale(", null, ",", null, ")");
      q.push({ i: i - 4, x: number_default(xa, xb) }, { i: i - 2, x: number_default(ya, yb) });
    } else if (xb !== 1 || yb !== 1) {
      s.push(pop(s) + "scale(" + xb + "," + yb + ")");
    }
  }
  return function(a, b) {
    var s = [], q = [];
    a = parse(a), b = parse(b);
    translate(a.translateX, a.translateY, b.translateX, b.translateY, s, q);
    rotate(a.rotate, b.rotate, s, q);
    skewX(a.skewX, b.skewX, s, q);
    scale(a.scaleX, a.scaleY, b.scaleX, b.scaleY, s, q);
    a = b = null;
    return function(t) {
      var i = -1, n = q.length, o;
      while (++i < n)
        s[(o = q[i]).i] = o.x(t);
      return s.join("");
    };
  };
}
var interpolateTransformCss = interpolateTransform(parseCss, "px, ", "px)", "deg)");
var interpolateTransformSvg = interpolateTransform(parseSvg, ", ", ")", ")");

// node_modules/d3-interpolate/src/zoom.js
var epsilon2 = 1e-12;
function cosh(x2) {
  return ((x2 = Math.exp(x2)) + 1 / x2) / 2;
}
function sinh(x2) {
  return ((x2 = Math.exp(x2)) - 1 / x2) / 2;
}
function tanh(x2) {
  return ((x2 = Math.exp(2 * x2)) - 1) / (x2 + 1);
}
var zoom_default = function zoomRho(rho, rho2, rho4) {
  function zoom(p0, p1) {
    var ux0 = p0[0], uy0 = p0[1], w0 = p0[2], ux1 = p1[0], uy1 = p1[1], w1 = p1[2], dx = ux1 - ux0, dy = uy1 - uy0, d2 = dx * dx + dy * dy, i, S;
    if (d2 < epsilon2) {
      S = Math.log(w1 / w0) / rho;
      i = function(t) {
        return [
          ux0 + t * dx,
          uy0 + t * dy,
          w0 * Math.exp(rho * t * S)
        ];
      };
    } else {
      var d1 = Math.sqrt(d2), b0 = (w1 * w1 - w0 * w0 + rho4 * d2) / (2 * w0 * rho2 * d1), b1 = (w1 * w1 - w0 * w0 - rho4 * d2) / (2 * w1 * rho2 * d1), r0 = Math.log(Math.sqrt(b0 * b0 + 1) - b0), r1 = Math.log(Math.sqrt(b1 * b1 + 1) - b1);
      S = (r1 - r0) / rho;
      i = function(t) {
        var s = t * S, coshr0 = cosh(r0), u = w0 / (rho2 * d1) * (coshr0 * tanh(rho * s + r0) - sinh(r0));
        return [
          ux0 + u * dx,
          uy0 + u * dy,
          w0 * coshr0 / cosh(rho * s + r0)
        ];
      };
    }
    i.duration = S * 1e3 * rho / Math.SQRT2;
    return i;
  }
  zoom.rho = function(_) {
    var _1 = Math.max(1e-3, +_), _2 = _1 * _1, _4 = _2 * _2;
    return zoomRho(_1, _2, _4);
  };
  return zoom;
}(Math.SQRT2, 2, 4);

// node_modules/d3-timer/src/timer.js
var frame = 0;
var timeout = 0;
var interval = 0;
var pokeDelay = 1e3;
var taskHead;
var taskTail;
var clockLast = 0;
var clockNow = 0;
var clockSkew = 0;
var clock = typeof performance === "object" && performance.now ? performance : Date;
var setFrame = typeof window === "object" && window.requestAnimationFrame ? window.requestAnimationFrame.bind(window) : function(f) {
  setTimeout(f, 17);
};
function now() {
  return clockNow || (setFrame(clearNow), clockNow = clock.now() + clockSkew);
}
function clearNow() {
  clockNow = 0;
}
function Timer() {
  this._call = this._time = this._next = null;
}
Timer.prototype = timer.prototype = {
  constructor: Timer,
  restart: function(callback, delay, time) {
    if (typeof callback !== "function")
      throw new TypeError("callback is not a function");
    time = (time == null ? now() : +time) + (delay == null ? 0 : +delay);
    if (!this._next && taskTail !== this) {
      if (taskTail)
        taskTail._next = this;
      else
        taskHead = this;
      taskTail = this;
    }
    this._call = callback;
    this._time = time;
    sleep();
  },
  stop: function() {
    if (this._call) {
      this._call = null;
      this._time = Infinity;
      sleep();
    }
  }
};
function timer(callback, delay, time) {
  var t = new Timer();
  t.restart(callback, delay, time);
  return t;
}
function timerFlush() {
  now();
  ++frame;
  var t = taskHead, e;
  while (t) {
    if ((e = clockNow - t._time) >= 0)
      t._call.call(void 0, e);
    t = t._next;
  }
  --frame;
}
function wake() {
  clockNow = (clockLast = clock.now()) + clockSkew;
  frame = timeout = 0;
  try {
    timerFlush();
  } finally {
    frame = 0;
    nap();
    clockNow = 0;
  }
}
function poke() {
  var now2 = clock.now(), delay = now2 - clockLast;
  if (delay > pokeDelay)
    clockSkew -= delay, clockLast = now2;
}
function nap() {
  var t0, t1 = taskHead, t2, time = Infinity;
  while (t1) {
    if (t1._call) {
      if (time > t1._time)
        time = t1._time;
      t0 = t1, t1 = t1._next;
    } else {
      t2 = t1._next, t1._next = null;
      t1 = t0 ? t0._next = t2 : taskHead = t2;
    }
  }
  taskTail = t0;
  sleep(time);
}
function sleep(time) {
  if (frame)
    return;
  if (timeout)
    timeout = clearTimeout(timeout);
  var delay = time - clockNow;
  if (delay > 24) {
    if (time < Infinity)
      timeout = setTimeout(wake, time - clock.now() - clockSkew);
    if (interval)
      interval = clearInterval(interval);
  } else {
    if (!interval)
      clockLast = clock.now(), interval = setInterval(poke, pokeDelay);
    frame = 1, setFrame(wake);
  }
}

// node_modules/d3-timer/src/timeout.js
function timeout_default(callback, delay, time) {
  var t = new Timer();
  delay = delay == null ? 0 : +delay;
  t.restart((elapsed) => {
    t.stop();
    callback(elapsed + delay);
  }, delay, time);
  return t;
}

// node_modules/d3-transition/src/transition/schedule.js
var emptyOn = dispatch_default("start", "end", "cancel", "interrupt");
var emptyTween = [];
var CREATED = 0;
var SCHEDULED = 1;
var STARTING = 2;
var STARTED = 3;
var RUNNING = 4;
var ENDING = 5;
var ENDED = 6;
function schedule_default(node, name, id2, index, group, timing) {
  var schedules = node.__transition;
  if (!schedules)
    node.__transition = {};
  else if (id2 in schedules)
    return;
  create(node, id2, {
    name,
    index,
    // For context during callback.
    group,
    // For context during callback.
    on: emptyOn,
    tween: emptyTween,
    time: timing.time,
    delay: timing.delay,
    duration: timing.duration,
    ease: timing.ease,
    timer: null,
    state: CREATED
  });
}
function init(node, id2) {
  var schedule = get2(node, id2);
  if (schedule.state > CREATED)
    throw new Error("too late; already scheduled");
  return schedule;
}
function set2(node, id2) {
  var schedule = get2(node, id2);
  if (schedule.state > STARTED)
    throw new Error("too late; already running");
  return schedule;
}
function get2(node, id2) {
  var schedule = node.__transition;
  if (!schedule || !(schedule = schedule[id2]))
    throw new Error("transition not found");
  return schedule;
}
function create(node, id2, self) {
  var schedules = node.__transition, tween;
  schedules[id2] = self;
  self.timer = timer(schedule, 0, self.time);
  function schedule(elapsed) {
    self.state = SCHEDULED;
    self.timer.restart(start2, self.delay, self.time);
    if (self.delay <= elapsed)
      start2(elapsed - self.delay);
  }
  function start2(elapsed) {
    var i, j, n, o;
    if (self.state !== SCHEDULED)
      return stop();
    for (i in schedules) {
      o = schedules[i];
      if (o.name !== self.name)
        continue;
      if (o.state === STARTED)
        return timeout_default(start2);
      if (o.state === RUNNING) {
        o.state = ENDED;
        o.timer.stop();
        o.on.call("interrupt", node, node.__data__, o.index, o.group);
        delete schedules[i];
      } else if (+i < id2) {
        o.state = ENDED;
        o.timer.stop();
        o.on.call("cancel", node, node.__data__, o.index, o.group);
        delete schedules[i];
      }
    }
    timeout_default(function() {
      if (self.state === STARTED) {
        self.state = RUNNING;
        self.timer.restart(tick, self.delay, self.time);
        tick(elapsed);
      }
    });
    self.state = STARTING;
    self.on.call("start", node, node.__data__, self.index, self.group);
    if (self.state !== STARTING)
      return;
    self.state = STARTED;
    tween = new Array(n = self.tween.length);
    for (i = 0, j = -1; i < n; ++i) {
      if (o = self.tween[i].value.call(node, node.__data__, self.index, self.group)) {
        tween[++j] = o;
      }
    }
    tween.length = j + 1;
  }
  function tick(elapsed) {
    var t = elapsed < self.duration ? self.ease.call(null, elapsed / self.duration) : (self.timer.restart(stop), self.state = ENDING, 1), i = -1, n = tween.length;
    while (++i < n) {
      tween[i].call(node, t);
    }
    if (self.state === ENDING) {
      self.on.call("end", node, node.__data__, self.index, self.group);
      stop();
    }
  }
  function stop() {
    self.state = ENDED;
    self.timer.stop();
    delete schedules[id2];
    for (var i in schedules)
      return;
    delete node.__transition;
  }
}

// node_modules/d3-transition/src/interrupt.js
function interrupt_default(node, name) {
  var schedules = node.__transition, schedule, active, empty2 = true, i;
  if (!schedules)
    return;
  name = name == null ? null : name + "";
  for (i in schedules) {
    if ((schedule = schedules[i]).name !== name) {
      empty2 = false;
      continue;
    }
    active = schedule.state > STARTING && schedule.state < ENDING;
    schedule.state = ENDED;
    schedule.timer.stop();
    schedule.on.call(active ? "interrupt" : "cancel", node, node.__data__, schedule.index, schedule.group);
    delete schedules[i];
  }
  if (empty2)
    delete node.__transition;
}

// node_modules/d3-transition/src/selection/interrupt.js
function interrupt_default2(name) {
  return this.each(function() {
    interrupt_default(this, name);
  });
}

// node_modules/d3-transition/src/transition/tween.js
function tweenRemove(id2, name) {
  var tween0, tween1;
  return function() {
    var schedule = set2(this, id2), tween = schedule.tween;
    if (tween !== tween0) {
      tween1 = tween0 = tween;
      for (var i = 0, n = tween1.length; i < n; ++i) {
        if (tween1[i].name === name) {
          tween1 = tween1.slice();
          tween1.splice(i, 1);
          break;
        }
      }
    }
    schedule.tween = tween1;
  };
}
function tweenFunction(id2, name, value) {
  var tween0, tween1;
  if (typeof value !== "function")
    throw new Error();
  return function() {
    var schedule = set2(this, id2), tween = schedule.tween;
    if (tween !== tween0) {
      tween1 = (tween0 = tween).slice();
      for (var t = { name, value }, i = 0, n = tween1.length; i < n; ++i) {
        if (tween1[i].name === name) {
          tween1[i] = t;
          break;
        }
      }
      if (i === n)
        tween1.push(t);
    }
    schedule.tween = tween1;
  };
}
function tween_default(name, value) {
  var id2 = this._id;
  name += "";
  if (arguments.length < 2) {
    var tween = get2(this.node(), id2).tween;
    for (var i = 0, n = tween.length, t; i < n; ++i) {
      if ((t = tween[i]).name === name) {
        return t.value;
      }
    }
    return null;
  }
  return this.each((value == null ? tweenRemove : tweenFunction)(id2, name, value));
}
function tweenValue(transition2, name, value) {
  var id2 = transition2._id;
  transition2.each(function() {
    var schedule = set2(this, id2);
    (schedule.value || (schedule.value = {}))[name] = value.apply(this, arguments);
  });
  return function(node) {
    return get2(node, id2).value[name];
  };
}

// node_modules/d3-transition/src/transition/interpolate.js
function interpolate_default(a, b) {
  var c;
  return (typeof b === "number" ? number_default : b instanceof color ? rgb_default : (c = color(b)) ? (b = c, rgb_default) : string_default)(a, b);
}

// node_modules/d3-transition/src/transition/attr.js
function attrRemove2(name) {
  return function() {
    this.removeAttribute(name);
  };
}
function attrRemoveNS2(fullname) {
  return function() {
    this.removeAttributeNS(fullname.space, fullname.local);
  };
}
function attrConstant2(name, interpolate, value1) {
  var string00, string1 = value1 + "", interpolate0;
  return function() {
    var string0 = this.getAttribute(name);
    return string0 === string1 ? null : string0 === string00 ? interpolate0 : interpolate0 = interpolate(string00 = string0, value1);
  };
}
function attrConstantNS2(fullname, interpolate, value1) {
  var string00, string1 = value1 + "", interpolate0;
  return function() {
    var string0 = this.getAttributeNS(fullname.space, fullname.local);
    return string0 === string1 ? null : string0 === string00 ? interpolate0 : interpolate0 = interpolate(string00 = string0, value1);
  };
}
function attrFunction2(name, interpolate, value) {
  var string00, string10, interpolate0;
  return function() {
    var string0, value1 = value(this), string1;
    if (value1 == null)
      return void this.removeAttribute(name);
    string0 = this.getAttribute(name);
    string1 = value1 + "";
    return string0 === string1 ? null : string0 === string00 && string1 === string10 ? interpolate0 : (string10 = string1, interpolate0 = interpolate(string00 = string0, value1));
  };
}
function attrFunctionNS2(fullname, interpolate, value) {
  var string00, string10, interpolate0;
  return function() {
    var string0, value1 = value(this), string1;
    if (value1 == null)
      return void this.removeAttributeNS(fullname.space, fullname.local);
    string0 = this.getAttributeNS(fullname.space, fullname.local);
    string1 = value1 + "";
    return string0 === string1 ? null : string0 === string00 && string1 === string10 ? interpolate0 : (string10 = string1, interpolate0 = interpolate(string00 = string0, value1));
  };
}
function attr_default2(name, value) {
  var fullname = namespace_default(name), i = fullname === "transform" ? interpolateTransformSvg : interpolate_default;
  return this.attrTween(name, typeof value === "function" ? (fullname.local ? attrFunctionNS2 : attrFunction2)(fullname, i, tweenValue(this, "attr." + name, value)) : value == null ? (fullname.local ? attrRemoveNS2 : attrRemove2)(fullname) : (fullname.local ? attrConstantNS2 : attrConstant2)(fullname, i, value));
}

// node_modules/d3-transition/src/transition/attrTween.js
function attrInterpolate(name, i) {
  return function(t) {
    this.setAttribute(name, i.call(this, t));
  };
}
function attrInterpolateNS(fullname, i) {
  return function(t) {
    this.setAttributeNS(fullname.space, fullname.local, i.call(this, t));
  };
}
function attrTweenNS(fullname, value) {
  var t0, i0;
  function tween() {
    var i = value.apply(this, arguments);
    if (i !== i0)
      t0 = (i0 = i) && attrInterpolateNS(fullname, i);
    return t0;
  }
  tween._value = value;
  return tween;
}
function attrTween(name, value) {
  var t0, i0;
  function tween() {
    var i = value.apply(this, arguments);
    if (i !== i0)
      t0 = (i0 = i) && attrInterpolate(name, i);
    return t0;
  }
  tween._value = value;
  return tween;
}
function attrTween_default(name, value) {
  var key = "attr." + name;
  if (arguments.length < 2)
    return (key = this.tween(key)) && key._value;
  if (value == null)
    return this.tween(key, null);
  if (typeof value !== "function")
    throw new Error();
  var fullname = namespace_default(name);
  return this.tween(key, (fullname.local ? attrTweenNS : attrTween)(fullname, value));
}

// node_modules/d3-transition/src/transition/delay.js
function delayFunction(id2, value) {
  return function() {
    init(this, id2).delay = +value.apply(this, arguments);
  };
}
function delayConstant(id2, value) {
  return value = +value, function() {
    init(this, id2).delay = value;
  };
}
function delay_default(value) {
  var id2 = this._id;
  return arguments.length ? this.each((typeof value === "function" ? delayFunction : delayConstant)(id2, value)) : get2(this.node(), id2).delay;
}

// node_modules/d3-transition/src/transition/duration.js
function durationFunction(id2, value) {
  return function() {
    set2(this, id2).duration = +value.apply(this, arguments);
  };
}
function durationConstant(id2, value) {
  return value = +value, function() {
    set2(this, id2).duration = value;
  };
}
function duration_default(value) {
  var id2 = this._id;
  return arguments.length ? this.each((typeof value === "function" ? durationFunction : durationConstant)(id2, value)) : get2(this.node(), id2).duration;
}

// node_modules/d3-transition/src/transition/ease.js
function easeConstant(id2, value) {
  if (typeof value !== "function")
    throw new Error();
  return function() {
    set2(this, id2).ease = value;
  };
}
function ease_default(value) {
  var id2 = this._id;
  return arguments.length ? this.each(easeConstant(id2, value)) : get2(this.node(), id2).ease;
}

// node_modules/d3-transition/src/transition/easeVarying.js
function easeVarying(id2, value) {
  return function() {
    var v = value.apply(this, arguments);
    if (typeof v !== "function")
      throw new Error();
    set2(this, id2).ease = v;
  };
}
function easeVarying_default(value) {
  if (typeof value !== "function")
    throw new Error();
  return this.each(easeVarying(this._id, value));
}

// node_modules/d3-transition/src/transition/filter.js
function filter_default2(match) {
  if (typeof match !== "function")
    match = matcher_default(match);
  for (var groups = this._groups, m = groups.length, subgroups = new Array(m), j = 0; j < m; ++j) {
    for (var group = groups[j], n = group.length, subgroup = subgroups[j] = [], node, i = 0; i < n; ++i) {
      if ((node = group[i]) && match.call(node, node.__data__, i, group)) {
        subgroup.push(node);
      }
    }
  }
  return new Transition(subgroups, this._parents, this._name, this._id);
}

// node_modules/d3-transition/src/transition/merge.js
function merge_default2(transition2) {
  if (transition2._id !== this._id)
    throw new Error();
  for (var groups0 = this._groups, groups1 = transition2._groups, m0 = groups0.length, m1 = groups1.length, m = Math.min(m0, m1), merges = new Array(m0), j = 0; j < m; ++j) {
    for (var group0 = groups0[j], group1 = groups1[j], n = group0.length, merge = merges[j] = new Array(n), node, i = 0; i < n; ++i) {
      if (node = group0[i] || group1[i]) {
        merge[i] = node;
      }
    }
  }
  for (; j < m0; ++j) {
    merges[j] = groups0[j];
  }
  return new Transition(merges, this._parents, this._name, this._id);
}

// node_modules/d3-transition/src/transition/on.js
function start(name) {
  return (name + "").trim().split(/^|\s+/).every(function(t) {
    var i = t.indexOf(".");
    if (i >= 0)
      t = t.slice(0, i);
    return !t || t === "start";
  });
}
function onFunction(id2, name, listener) {
  var on0, on1, sit = start(name) ? init : set2;
  return function() {
    var schedule = sit(this, id2), on = schedule.on;
    if (on !== on0)
      (on1 = (on0 = on).copy()).on(name, listener);
    schedule.on = on1;
  };
}
function on_default2(name, listener) {
  var id2 = this._id;
  return arguments.length < 2 ? get2(this.node(), id2).on.on(name) : this.each(onFunction(id2, name, listener));
}

// node_modules/d3-transition/src/transition/remove.js
function removeFunction(id2) {
  return function() {
    var parent = this.parentNode;
    for (var i in this.__transition)
      if (+i !== id2)
        return;
    if (parent)
      parent.removeChild(this);
  };
}
function remove_default2() {
  return this.on("end.remove", removeFunction(this._id));
}

// node_modules/d3-transition/src/transition/select.js
function select_default3(select) {
  var name = this._name, id2 = this._id;
  if (typeof select !== "function")
    select = selector_default(select);
  for (var groups = this._groups, m = groups.length, subgroups = new Array(m), j = 0; j < m; ++j) {
    for (var group = groups[j], n = group.length, subgroup = subgroups[j] = new Array(n), node, subnode, i = 0; i < n; ++i) {
      if ((node = group[i]) && (subnode = select.call(node, node.__data__, i, group))) {
        if ("__data__" in node)
          subnode.__data__ = node.__data__;
        subgroup[i] = subnode;
        schedule_default(subgroup[i], name, id2, i, subgroup, get2(node, id2));
      }
    }
  }
  return new Transition(subgroups, this._parents, name, id2);
}

// node_modules/d3-transition/src/transition/selectAll.js
function selectAll_default2(select) {
  var name = this._name, id2 = this._id;
  if (typeof select !== "function")
    select = selectorAll_default(select);
  for (var groups = this._groups, m = groups.length, subgroups = [], parents = [], j = 0; j < m; ++j) {
    for (var group = groups[j], n = group.length, node, i = 0; i < n; ++i) {
      if (node = group[i]) {
        for (var children2 = select.call(node, node.__data__, i, group), child, inherit2 = get2(node, id2), k = 0, l = children2.length; k < l; ++k) {
          if (child = children2[k]) {
            schedule_default(child, name, id2, k, children2, inherit2);
          }
        }
        subgroups.push(children2);
        parents.push(node);
      }
    }
  }
  return new Transition(subgroups, parents, name, id2);
}

// node_modules/d3-transition/src/transition/selection.js
var Selection2 = selection_default.prototype.constructor;
function selection_default2() {
  return new Selection2(this._groups, this._parents);
}

// node_modules/d3-transition/src/transition/style.js
function styleNull(name, interpolate) {
  var string00, string10, interpolate0;
  return function() {
    var string0 = styleValue(this, name), string1 = (this.style.removeProperty(name), styleValue(this, name));
    return string0 === string1 ? null : string0 === string00 && string1 === string10 ? interpolate0 : interpolate0 = interpolate(string00 = string0, string10 = string1);
  };
}
function styleRemove2(name) {
  return function() {
    this.style.removeProperty(name);
  };
}
function styleConstant2(name, interpolate, value1) {
  var string00, string1 = value1 + "", interpolate0;
  return function() {
    var string0 = styleValue(this, name);
    return string0 === string1 ? null : string0 === string00 ? interpolate0 : interpolate0 = interpolate(string00 = string0, value1);
  };
}
function styleFunction2(name, interpolate, value) {
  var string00, string10, interpolate0;
  return function() {
    var string0 = styleValue(this, name), value1 = value(this), string1 = value1 + "";
    if (value1 == null)
      string1 = value1 = (this.style.removeProperty(name), styleValue(this, name));
    return string0 === string1 ? null : string0 === string00 && string1 === string10 ? interpolate0 : (string10 = string1, interpolate0 = interpolate(string00 = string0, value1));
  };
}
function styleMaybeRemove(id2, name) {
  var on0, on1, listener0, key = "style." + name, event = "end." + key, remove2;
  return function() {
    var schedule = set2(this, id2), on = schedule.on, listener = schedule.value[key] == null ? remove2 || (remove2 = styleRemove2(name)) : void 0;
    if (on !== on0 || listener0 !== listener)
      (on1 = (on0 = on).copy()).on(event, listener0 = listener);
    schedule.on = on1;
  };
}
function style_default2(name, value, priority) {
  var i = (name += "") === "transform" ? interpolateTransformCss : interpolate_default;
  return value == null ? this.styleTween(name, styleNull(name, i)).on("end.style." + name, styleRemove2(name)) : typeof value === "function" ? this.styleTween(name, styleFunction2(name, i, tweenValue(this, "style." + name, value))).each(styleMaybeRemove(this._id, name)) : this.styleTween(name, styleConstant2(name, i, value), priority).on("end.style." + name, null);
}

// node_modules/d3-transition/src/transition/styleTween.js
function styleInterpolate(name, i, priority) {
  return function(t) {
    this.style.setProperty(name, i.call(this, t), priority);
  };
}
function styleTween(name, value, priority) {
  var t, i0;
  function tween() {
    var i = value.apply(this, arguments);
    if (i !== i0)
      t = (i0 = i) && styleInterpolate(name, i, priority);
    return t;
  }
  tween._value = value;
  return tween;
}
function styleTween_default(name, value, priority) {
  var key = "style." + (name += "");
  if (arguments.length < 2)
    return (key = this.tween(key)) && key._value;
  if (value == null)
    return this.tween(key, null);
  if (typeof value !== "function")
    throw new Error();
  return this.tween(key, styleTween(name, value, priority == null ? "" : priority));
}

// node_modules/d3-transition/src/transition/text.js
function textConstant2(value) {
  return function() {
    this.textContent = value;
  };
}
function textFunction2(value) {
  return function() {
    var value1 = value(this);
    this.textContent = value1 == null ? "" : value1;
  };
}
function text_default2(value) {
  return this.tween("text", typeof value === "function" ? textFunction2(tweenValue(this, "text", value)) : textConstant2(value == null ? "" : value + ""));
}

// node_modules/d3-transition/src/transition/textTween.js
function textInterpolate(i) {
  return function(t) {
    this.textContent = i.call(this, t);
  };
}
function textTween(value) {
  var t0, i0;
  function tween() {
    var i = value.apply(this, arguments);
    if (i !== i0)
      t0 = (i0 = i) && textInterpolate(i);
    return t0;
  }
  tween._value = value;
  return tween;
}
function textTween_default(value) {
  var key = "text";
  if (arguments.length < 1)
    return (key = this.tween(key)) && key._value;
  if (value == null)
    return this.tween(key, null);
  if (typeof value !== "function")
    throw new Error();
  return this.tween(key, textTween(value));
}

// node_modules/d3-transition/src/transition/transition.js
function transition_default() {
  var name = this._name, id0 = this._id, id1 = newId();
  for (var groups = this._groups, m = groups.length, j = 0; j < m; ++j) {
    for (var group = groups[j], n = group.length, node, i = 0; i < n; ++i) {
      if (node = group[i]) {
        var inherit2 = get2(node, id0);
        schedule_default(node, name, id1, i, group, {
          time: inherit2.time + inherit2.delay + inherit2.duration,
          delay: 0,
          duration: inherit2.duration,
          ease: inherit2.ease
        });
      }
    }
  }
  return new Transition(groups, this._parents, name, id1);
}

// node_modules/d3-transition/src/transition/end.js
function end_default() {
  var on0, on1, that = this, id2 = that._id, size = that.size();
  return new Promise(function(resolve, reject) {
    var cancel = { value: reject }, end = { value: function() {
      if (--size === 0)
        resolve();
    } };
    that.each(function() {
      var schedule = set2(this, id2), on = schedule.on;
      if (on !== on0) {
        on1 = (on0 = on).copy();
        on1._.cancel.push(cancel);
        on1._.interrupt.push(cancel);
        on1._.end.push(end);
      }
      schedule.on = on1;
    });
    if (size === 0)
      resolve();
  });
}

// node_modules/d3-transition/src/transition/index.js
var id = 0;
function Transition(groups, parents, name, id2) {
  this._groups = groups;
  this._parents = parents;
  this._name = name;
  this._id = id2;
}
function transition(name) {
  return selection_default().transition(name);
}
function newId() {
  return ++id;
}
var selection_prototype = selection_default.prototype;
Transition.prototype = transition.prototype = {
  constructor: Transition,
  select: select_default3,
  selectAll: selectAll_default2,
  selectChild: selection_prototype.selectChild,
  selectChildren: selection_prototype.selectChildren,
  filter: filter_default2,
  merge: merge_default2,
  selection: selection_default2,
  transition: transition_default,
  call: selection_prototype.call,
  nodes: selection_prototype.nodes,
  node: selection_prototype.node,
  size: selection_prototype.size,
  empty: selection_prototype.empty,
  each: selection_prototype.each,
  on: on_default2,
  attr: attr_default2,
  attrTween: attrTween_default,
  style: style_default2,
  styleTween: styleTween_default,
  text: text_default2,
  textTween: textTween_default,
  remove: remove_default2,
  tween: tween_default,
  delay: delay_default,
  duration: duration_default,
  ease: ease_default,
  easeVarying: easeVarying_default,
  end: end_default,
  [Symbol.iterator]: selection_prototype[Symbol.iterator]
};

// node_modules/d3-ease/src/cubic.js
function cubicInOut(t) {
  return ((t *= 2) <= 1 ? t * t * t : (t -= 2) * t * t + 2) / 2;
}

// node_modules/d3-transition/src/selection/transition.js
var defaultTiming = {
  time: null,
  // Set on use.
  delay: 0,
  duration: 250,
  ease: cubicInOut
};
function inherit(node, id2) {
  var timing;
  while (!(timing = node.__transition) || !(timing = timing[id2])) {
    if (!(node = node.parentNode)) {
      throw new Error(`transition ${id2} not found`);
    }
  }
  return timing;
}
function transition_default2(name) {
  var id2, timing;
  if (name instanceof Transition) {
    id2 = name._id, name = name._name;
  } else {
    id2 = newId(), (timing = defaultTiming).time = now(), name = name == null ? null : name + "";
  }
  for (var groups = this._groups, m = groups.length, j = 0; j < m; ++j) {
    for (var group = groups[j], n = group.length, node, i = 0; i < n; ++i) {
      if (node = group[i]) {
        schedule_default(node, name, id2, i, group, timing || inherit(node, id2));
      }
    }
  }
  return new Transition(groups, this._parents, name, id2);
}

// node_modules/d3-transition/src/selection/index.js
selection_default.prototype.interrupt = interrupt_default2;
selection_default.prototype.transition = transition_default2;

// node_modules/d3-brush/src/brush.js
var { abs, max, min } = Math;
function number1(e) {
  return [+e[0], +e[1]];
}
function number2(e) {
  return [number1(e[0]), number1(e[1])];
}
var X = {
  name: "x",
  handles: ["w", "e"].map(type),
  input: function(x2, e) {
    return x2 == null ? null : [[+x2[0], e[0][1]], [+x2[1], e[1][1]]];
  },
  output: function(xy) {
    return xy && [xy[0][0], xy[1][0]];
  }
};
var Y = {
  name: "y",
  handles: ["n", "s"].map(type),
  input: function(y2, e) {
    return y2 == null ? null : [[e[0][0], +y2[0]], [e[1][0], +y2[1]]];
  },
  output: function(xy) {
    return xy && [xy[0][1], xy[1][1]];
  }
};
var XY = {
  name: "xy",
  handles: ["n", "w", "e", "s", "nw", "ne", "sw", "se"].map(type),
  input: function(xy) {
    return xy == null ? null : number2(xy);
  },
  output: function(xy) {
    return xy;
  }
};
function type(t) {
  return { type: t };
}

// node_modules/d3-path/src/path.js
var pi = Math.PI;
var tau = 2 * pi;
var epsilon = 1e-6;
var tauEpsilon = tau - epsilon;
function append(strings) {
  this._ += strings[0];
  for (let i = 1, n = strings.length; i < n; ++i) {
    this._ += arguments[i] + strings[i];
  }
}
function appendRound(digits) {
  let d = Math.floor(digits);
  if (!(d >= 0))
    throw new Error(`invalid digits: ${digits}`);
  if (d > 15)
    return append;
  const k = 10 ** d;
  return function(strings) {
    this._ += strings[0];
    for (let i = 1, n = strings.length; i < n; ++i) {
      this._ += Math.round(arguments[i] * k) / k + strings[i];
    }
  };
}
var Path = class {
  constructor(digits) {
    this._x0 = this._y0 = // start of current subpath
    this._x1 = this._y1 = null;
    this._ = "";
    this._append = digits == null ? append : appendRound(digits);
  }
  moveTo(x2, y2) {
    this._append`M${this._x0 = this._x1 = +x2},${this._y0 = this._y1 = +y2}`;
  }
  closePath() {
    if (this._x1 !== null) {
      this._x1 = this._x0, this._y1 = this._y0;
      this._append`Z`;
    }
  }
  lineTo(x2, y2) {
    this._append`L${this._x1 = +x2},${this._y1 = +y2}`;
  }
  quadraticCurveTo(x1, y1, x2, y2) {
    this._append`Q${+x1},${+y1},${this._x1 = +x2},${this._y1 = +y2}`;
  }
  bezierCurveTo(x1, y1, x2, y2, x3, y3) {
    this._append`C${+x1},${+y1},${+x2},${+y2},${this._x1 = +x3},${this._y1 = +y3}`;
  }
  arcTo(x1, y1, x2, y2, r) {
    x1 = +x1, y1 = +y1, x2 = +x2, y2 = +y2, r = +r;
    if (r < 0)
      throw new Error(`negative radius: ${r}`);
    let x0 = this._x1, y0 = this._y1, x21 = x2 - x1, y21 = y2 - y1, x01 = x0 - x1, y01 = y0 - y1, l01_2 = x01 * x01 + y01 * y01;
    if (this._x1 === null) {
      this._append`M${this._x1 = x1},${this._y1 = y1}`;
    } else if (!(l01_2 > epsilon))
      ;
    else if (!(Math.abs(y01 * x21 - y21 * x01) > epsilon) || !r) {
      this._append`L${this._x1 = x1},${this._y1 = y1}`;
    } else {
      let x20 = x2 - x0, y20 = y2 - y0, l21_2 = x21 * x21 + y21 * y21, l20_2 = x20 * x20 + y20 * y20, l21 = Math.sqrt(l21_2), l01 = Math.sqrt(l01_2), l = r * Math.tan((pi - Math.acos((l21_2 + l01_2 - l20_2) / (2 * l21 * l01))) / 2), t01 = l / l01, t21 = l / l21;
      if (Math.abs(t01 - 1) > epsilon) {
        this._append`L${x1 + t01 * x01},${y1 + t01 * y01}`;
      }
      this._append`A${r},${r},0,0,${+(y01 * x20 > x01 * y20)},${this._x1 = x1 + t21 * x21},${this._y1 = y1 + t21 * y21}`;
    }
  }
  arc(x2, y2, r, a0, a1, ccw) {
    x2 = +x2, y2 = +y2, r = +r, ccw = !!ccw;
    if (r < 0)
      throw new Error(`negative radius: ${r}`);
    let dx = r * Math.cos(a0), dy = r * Math.sin(a0), x0 = x2 + dx, y0 = y2 + dy, cw = 1 ^ ccw, da = ccw ? a0 - a1 : a1 - a0;
    if (this._x1 === null) {
      this._append`M${x0},${y0}`;
    } else if (Math.abs(this._x1 - x0) > epsilon || Math.abs(this._y1 - y0) > epsilon) {
      this._append`L${x0},${y0}`;
    }
    if (!r)
      return;
    if (da < 0)
      da = da % tau + tau;
    if (da > tauEpsilon) {
      this._append`A${r},${r},0,1,${cw},${x2 - dx},${y2 - dy}A${r},${r},0,1,${cw},${this._x1 = x0},${this._y1 = y0}`;
    } else if (da > epsilon) {
      this._append`A${r},${r},0,${+(da >= pi)},${cw},${this._x1 = x2 + r * Math.cos(a1)},${this._y1 = y2 + r * Math.sin(a1)}`;
    }
  }
  rect(x2, y2, w, h) {
    this._append`M${this._x0 = this._x1 = +x2},${this._y0 = this._y1 = +y2}h${w = +w}v${+h}h${-w}Z`;
  }
  toString() {
    return this._;
  }
};
function path() {
  return new Path();
}
path.prototype = Path.prototype;

// node_modules/d3-format/src/formatDecimal.js
function formatDecimal_default(x2) {
  return Math.abs(x2 = Math.round(x2)) >= 1e21 ? x2.toLocaleString("en").replace(/,/g, "") : x2.toString(10);
}
function formatDecimalParts(x2, p) {
  if (!isFinite(x2) || x2 === 0)
    return null;
  var i = (x2 = p ? x2.toExponential(p - 1) : x2.toExponential()).indexOf("e"), coefficient = x2.slice(0, i);
  return [
    coefficient.length > 1 ? coefficient[0] + coefficient.slice(2) : coefficient,
    +x2.slice(i + 1)
  ];
}

// node_modules/d3-format/src/exponent.js
function exponent_default(x2) {
  return x2 = formatDecimalParts(Math.abs(x2)), x2 ? x2[1] : NaN;
}

// node_modules/d3-format/src/formatGroup.js
function formatGroup_default(grouping, thousands) {
  return function(value, width) {
    var i = value.length, t = [], j = 0, g = grouping[0], length = 0;
    while (i > 0 && g > 0) {
      if (length + g + 1 > width)
        g = Math.max(1, width - length);
      t.push(value.substring(i -= g, i + g));
      if ((length += g + 1) > width)
        break;
      g = grouping[j = (j + 1) % grouping.length];
    }
    return t.reverse().join(thousands);
  };
}

// node_modules/d3-format/src/formatNumerals.js
function formatNumerals_default(numerals) {
  return function(value) {
    return value.replace(/[0-9]/g, function(i) {
      return numerals[+i];
    });
  };
}

// node_modules/d3-format/src/formatSpecifier.js
var re = /^(?:(.)?([<>=^]))?([+\-( ])?([$#])?(0)?(\d+)?(,)?(\.\d+)?(~)?([a-z%])?$/i;
function formatSpecifier(specifier) {
  if (!(match = re.exec(specifier)))
    throw new Error("invalid format: " + specifier);
  var match;
  return new FormatSpecifier({
    fill: match[1],
    align: match[2],
    sign: match[3],
    symbol: match[4],
    zero: match[5],
    width: match[6],
    comma: match[7],
    precision: match[8] && match[8].slice(1),
    trim: match[9],
    type: match[10]
  });
}
formatSpecifier.prototype = FormatSpecifier.prototype;
function FormatSpecifier(specifier) {
  this.fill = specifier.fill === void 0 ? " " : specifier.fill + "";
  this.align = specifier.align === void 0 ? ">" : specifier.align + "";
  this.sign = specifier.sign === void 0 ? "-" : specifier.sign + "";
  this.symbol = specifier.symbol === void 0 ? "" : specifier.symbol + "";
  this.zero = !!specifier.zero;
  this.width = specifier.width === void 0 ? void 0 : +specifier.width;
  this.comma = !!specifier.comma;
  this.precision = specifier.precision === void 0 ? void 0 : +specifier.precision;
  this.trim = !!specifier.trim;
  this.type = specifier.type === void 0 ? "" : specifier.type + "";
}
FormatSpecifier.prototype.toString = function() {
  return this.fill + this.align + this.sign + this.symbol + (this.zero ? "0" : "") + (this.width === void 0 ? "" : Math.max(1, this.width | 0)) + (this.comma ? "," : "") + (this.precision === void 0 ? "" : "." + Math.max(0, this.precision | 0)) + (this.trim ? "~" : "") + this.type;
};

// node_modules/d3-format/src/formatTrim.js
function formatTrim_default(s) {
  out:
    for (var n = s.length, i = 1, i0 = -1, i1; i < n; ++i) {
      switch (s[i]) {
        case ".":
          i0 = i1 = i;
          break;
        case "0":
          if (i0 === 0)
            i0 = i;
          i1 = i;
          break;
        default:
          if (!+s[i])
            break out;
          if (i0 > 0)
            i0 = 0;
          break;
      }
    }
  return i0 > 0 ? s.slice(0, i0) + s.slice(i1 + 1) : s;
}

// node_modules/d3-format/src/formatPrefixAuto.js
var prefixExponent;
function formatPrefixAuto_default(x2, p) {
  var d = formatDecimalParts(x2, p);
  if (!d)
    return prefixExponent = void 0, x2.toPrecision(p);
  var coefficient = d[0], exponent = d[1], i = exponent - (prefixExponent = Math.max(-8, Math.min(8, Math.floor(exponent / 3))) * 3) + 1, n = coefficient.length;
  return i === n ? coefficient : i > n ? coefficient + new Array(i - n + 1).join("0") : i > 0 ? coefficient.slice(0, i) + "." + coefficient.slice(i) : "0." + new Array(1 - i).join("0") + formatDecimalParts(x2, Math.max(0, p + i - 1))[0];
}

// node_modules/d3-format/src/formatRounded.js
function formatRounded_default(x2, p) {
  var d = formatDecimalParts(x2, p);
  if (!d)
    return x2 + "";
  var coefficient = d[0], exponent = d[1];
  return exponent < 0 ? "0." + new Array(-exponent).join("0") + coefficient : coefficient.length > exponent + 1 ? coefficient.slice(0, exponent + 1) + "." + coefficient.slice(exponent + 1) : coefficient + new Array(exponent - coefficient.length + 2).join("0");
}

// node_modules/d3-format/src/formatTypes.js
var formatTypes_default = {
  "%": (x2, p) => (x2 * 100).toFixed(p),
  "b": (x2) => Math.round(x2).toString(2),
  "c": (x2) => x2 + "",
  "d": formatDecimal_default,
  "e": (x2, p) => x2.toExponential(p),
  "f": (x2, p) => x2.toFixed(p),
  "g": (x2, p) => x2.toPrecision(p),
  "o": (x2) => Math.round(x2).toString(8),
  "p": (x2, p) => formatRounded_default(x2 * 100, p),
  "r": formatRounded_default,
  "s": formatPrefixAuto_default,
  "X": (x2) => Math.round(x2).toString(16).toUpperCase(),
  "x": (x2) => Math.round(x2).toString(16)
};

// node_modules/d3-format/src/identity.js
function identity_default(x2) {
  return x2;
}

// node_modules/d3-format/src/locale.js
var map = Array.prototype.map;
var prefixes = ["y", "z", "a", "f", "p", "n", "\xB5", "m", "", "k", "M", "G", "T", "P", "E", "Z", "Y"];
function locale_default(locale2) {
  var group = locale2.grouping === void 0 || locale2.thousands === void 0 ? identity_default : formatGroup_default(map.call(locale2.grouping, Number), locale2.thousands + ""), currencyPrefix = locale2.currency === void 0 ? "" : locale2.currency[0] + "", currencySuffix = locale2.currency === void 0 ? "" : locale2.currency[1] + "", decimal = locale2.decimal === void 0 ? "." : locale2.decimal + "", numerals = locale2.numerals === void 0 ? identity_default : formatNumerals_default(map.call(locale2.numerals, String)), percent = locale2.percent === void 0 ? "%" : locale2.percent + "", minus = locale2.minus === void 0 ? "\u2212" : locale2.minus + "", nan = locale2.nan === void 0 ? "NaN" : locale2.nan + "";
  function newFormat(specifier, options) {
    specifier = formatSpecifier(specifier);
    var fill = specifier.fill, align = specifier.align, sign = specifier.sign, symbol = specifier.symbol, zero3 = specifier.zero, width = specifier.width, comma = specifier.comma, precision = specifier.precision, trim = specifier.trim, type2 = specifier.type;
    if (type2 === "n")
      comma = true, type2 = "g";
    else if (!formatTypes_default[type2])
      precision === void 0 && (precision = 12), trim = true, type2 = "g";
    if (zero3 || fill === "0" && align === "=")
      zero3 = true, fill = "0", align = "=";
    var prefix = (options && options.prefix !== void 0 ? options.prefix : "") + (symbol === "$" ? currencyPrefix : symbol === "#" && /[boxX]/.test(type2) ? "0" + type2.toLowerCase() : ""), suffix = (symbol === "$" ? currencySuffix : /[%p]/.test(type2) ? percent : "") + (options && options.suffix !== void 0 ? options.suffix : "");
    var formatType = formatTypes_default[type2], maybeSuffix = /[defgprs%]/.test(type2);
    precision = precision === void 0 ? 6 : /[gprs]/.test(type2) ? Math.max(1, Math.min(21, precision)) : Math.max(0, Math.min(20, precision));
    function format2(value) {
      var valuePrefix = prefix, valueSuffix = suffix, i, n, c;
      if (type2 === "c") {
        valueSuffix = formatType(value) + valueSuffix;
        value = "";
      } else {
        value = +value;
        var valueNegative = value < 0 || 1 / value < 0;
        value = isNaN(value) ? nan : formatType(Math.abs(value), precision);
        if (trim)
          value = formatTrim_default(value);
        if (valueNegative && +value === 0 && sign !== "+")
          valueNegative = false;
        valuePrefix = (valueNegative ? sign === "(" ? sign : minus : sign === "-" || sign === "(" ? "" : sign) + valuePrefix;
        valueSuffix = (type2 === "s" && !isNaN(value) && prefixExponent !== void 0 ? prefixes[8 + prefixExponent / 3] : "") + valueSuffix + (valueNegative && sign === "(" ? ")" : "");
        if (maybeSuffix) {
          i = -1, n = value.length;
          while (++i < n) {
            if (c = value.charCodeAt(i), 48 > c || c > 57) {
              valueSuffix = (c === 46 ? decimal + value.slice(i + 1) : value.slice(i)) + valueSuffix;
              value = value.slice(0, i);
              break;
            }
          }
        }
      }
      if (comma && !zero3)
        value = group(value, Infinity);
      var length = valuePrefix.length + value.length + valueSuffix.length, padding = length < width ? new Array(width - length + 1).join(fill) : "";
      if (comma && zero3)
        value = group(padding + value, padding.length ? width - valueSuffix.length : Infinity), padding = "";
      switch (align) {
        case "<":
          value = valuePrefix + value + valueSuffix + padding;
          break;
        case "=":
          value = valuePrefix + padding + value + valueSuffix;
          break;
        case "^":
          value = padding.slice(0, length = padding.length >> 1) + valuePrefix + value + valueSuffix + padding.slice(length);
          break;
        default:
          value = padding + valuePrefix + value + valueSuffix;
          break;
      }
      return numerals(value);
    }
    format2.toString = function() {
      return specifier + "";
    };
    return format2;
  }
  function formatPrefix2(specifier, value) {
    var e = Math.max(-8, Math.min(8, Math.floor(exponent_default(value) / 3))) * 3, k = Math.pow(10, -e), f = newFormat((specifier = formatSpecifier(specifier), specifier.type = "f", specifier), { suffix: prefixes[8 + e / 3] });
    return function(value2) {
      return f(k * value2);
    };
  }
  return {
    format: newFormat,
    formatPrefix: formatPrefix2
  };
}

// node_modules/d3-format/src/defaultLocale.js
var locale;
var format;
var formatPrefix;
defaultLocale({
  thousands: ",",
  grouping: [3],
  currency: ["$", ""]
});
function defaultLocale(definition) {
  locale = locale_default(definition);
  format = locale.format;
  formatPrefix = locale.formatPrefix;
  return locale;
}

// node_modules/d3-format/src/precisionFixed.js
function precisionFixed_default(step) {
  return Math.max(0, -exponent_default(Math.abs(step)));
}

// node_modules/d3-format/src/precisionPrefix.js
function precisionPrefix_default(step, value) {
  return Math.max(0, Math.max(-8, Math.min(8, Math.floor(exponent_default(value) / 3))) * 3 - exponent_default(Math.abs(step)));
}

// node_modules/d3-format/src/precisionRound.js
function precisionRound_default(step, max2) {
  step = Math.abs(step), max2 = Math.abs(max2) - step;
  return Math.max(0, exponent_default(max2) - exponent_default(step)) + 1;
}

// node_modules/d3-scale/src/init.js
function initRange(domain, range) {
  switch (arguments.length) {
    case 0:
      break;
    case 1:
      this.range(domain);
      break;
    default:
      this.range(range).domain(domain);
      break;
  }
  return this;
}

// node_modules/d3-scale/src/constant.js
function constants(x2) {
  return function() {
    return x2;
  };
}

// node_modules/d3-scale/src/number.js
function number3(x2) {
  return +x2;
}

// node_modules/d3-scale/src/continuous.js
var unit = [0, 1];
function identity2(x2) {
  return x2;
}
function normalize(a, b) {
  return (b -= a = +a) ? function(x2) {
    return (x2 - a) / b;
  } : constants(isNaN(b) ? NaN : 0.5);
}
function clamper(a, b) {
  var t;
  if (a > b)
    t = a, a = b, b = t;
  return function(x2) {
    return Math.max(a, Math.min(b, x2));
  };
}
function bimap(domain, range, interpolate) {
  var d0 = domain[0], d1 = domain[1], r0 = range[0], r1 = range[1];
  if (d1 < d0)
    d0 = normalize(d1, d0), r0 = interpolate(r1, r0);
  else
    d0 = normalize(d0, d1), r0 = interpolate(r0, r1);
  return function(x2) {
    return r0(d0(x2));
  };
}
function polymap(domain, range, interpolate) {
  var j = Math.min(domain.length, range.length) - 1, d = new Array(j), r = new Array(j), i = -1;
  if (domain[j] < domain[0]) {
    domain = domain.slice().reverse();
    range = range.slice().reverse();
  }
  while (++i < j) {
    d[i] = normalize(domain[i], domain[i + 1]);
    r[i] = interpolate(range[i], range[i + 1]);
  }
  return function(x2) {
    var i2 = bisect_default(domain, x2, 1, j) - 1;
    return r[i2](d[i2](x2));
  };
}
function copy(source, target) {
  return target.domain(source.domain()).range(source.range()).interpolate(source.interpolate()).clamp(source.clamp()).unknown(source.unknown());
}
function transformer() {
  var domain = unit, range = unit, interpolate = value_default, transform2, untransform, unknown, clamp = identity2, piecewise, output, input;
  function rescale() {
    var n = Math.min(domain.length, range.length);
    if (clamp !== identity2)
      clamp = clamper(domain[0], domain[n - 1]);
    piecewise = n > 2 ? polymap : bimap;
    output = input = null;
    return scale;
  }
  function scale(x2) {
    return x2 == null || isNaN(x2 = +x2) ? unknown : (output || (output = piecewise(domain.map(transform2), range, interpolate)))(transform2(clamp(x2)));
  }
  scale.invert = function(y2) {
    return clamp(untransform((input || (input = piecewise(range, domain.map(transform2), number_default)))(y2)));
  };
  scale.domain = function(_) {
    return arguments.length ? (domain = Array.from(_, number3), rescale()) : domain.slice();
  };
  scale.range = function(_) {
    return arguments.length ? (range = Array.from(_), rescale()) : range.slice();
  };
  scale.rangeRound = function(_) {
    return range = Array.from(_), interpolate = round_default, rescale();
  };
  scale.clamp = function(_) {
    return arguments.length ? (clamp = _ ? true : identity2, rescale()) : clamp !== identity2;
  };
  scale.interpolate = function(_) {
    return arguments.length ? (interpolate = _, rescale()) : interpolate;
  };
  scale.unknown = function(_) {
    return arguments.length ? (unknown = _, scale) : unknown;
  };
  return function(t, u) {
    transform2 = t, untransform = u;
    return rescale();
  };
}
function continuous() {
  return transformer()(identity2, identity2);
}

// node_modules/d3-scale/src/tickFormat.js
function tickFormat(start2, stop, count, specifier) {
  var step = tickStep(start2, stop, count), precision;
  specifier = formatSpecifier(specifier == null ? ",f" : specifier);
  switch (specifier.type) {
    case "s": {
      var value = Math.max(Math.abs(start2), Math.abs(stop));
      if (specifier.precision == null && !isNaN(precision = precisionPrefix_default(step, value)))
        specifier.precision = precision;
      return formatPrefix(specifier, value);
    }
    case "":
    case "e":
    case "g":
    case "p":
    case "r": {
      if (specifier.precision == null && !isNaN(precision = precisionRound_default(step, Math.max(Math.abs(start2), Math.abs(stop)))))
        specifier.precision = precision - (specifier.type === "e");
      break;
    }
    case "f":
    case "%": {
      if (specifier.precision == null && !isNaN(precision = precisionFixed_default(step)))
        specifier.precision = precision - (specifier.type === "%") * 2;
      break;
    }
  }
  return format(specifier);
}

// node_modules/d3-scale/src/linear.js
function linearish(scale) {
  var domain = scale.domain;
  scale.ticks = function(count) {
    var d = domain();
    return ticks(d[0], d[d.length - 1], count == null ? 10 : count);
  };
  scale.tickFormat = function(count, specifier) {
    var d = domain();
    return tickFormat(d[0], d[d.length - 1], count == null ? 10 : count, specifier);
  };
  scale.nice = function(count) {
    if (count == null)
      count = 10;
    var d = domain();
    var i0 = 0;
    var i1 = d.length - 1;
    var start2 = d[i0];
    var stop = d[i1];
    var prestep;
    var step;
    var maxIter = 10;
    if (stop < start2) {
      step = start2, start2 = stop, stop = step;
      step = i0, i0 = i1, i1 = step;
    }
    while (maxIter-- > 0) {
      step = tickIncrement(start2, stop, count);
      if (step === prestep) {
        d[i0] = start2;
        d[i1] = stop;
        return domain(d);
      } else if (step > 0) {
        start2 = Math.floor(start2 / step) * step;
        stop = Math.ceil(stop / step) * step;
      } else if (step < 0) {
        start2 = Math.ceil(start2 * step) / step;
        stop = Math.floor(stop * step) / step;
      } else {
        break;
      }
      prestep = step;
    }
    return scale;
  };
  return scale;
}
function linear2() {
  var scale = continuous();
  scale.copy = function() {
    return copy(scale, linear2());
  };
  initRange.apply(scale, arguments);
  return linearish(scale);
}

// node_modules/d3-shape/src/constant.js
function constant_default4(x2) {
  return function constant() {
    return x2;
  };
}

// node_modules/d3-shape/src/math.js
var sqrt = Math.sqrt;
var epsilon3 = 1e-12;
var pi2 = Math.PI;
var halfPi = pi2 / 2;
var tau2 = 2 * pi2;

// node_modules/d3-shape/src/path.js
function withPath(shape) {
  let digits = 3;
  shape.digits = function(_) {
    if (!arguments.length)
      return digits;
    if (_ == null) {
      digits = null;
    } else {
      const d = Math.floor(_);
      if (!(d >= 0))
        throw new RangeError(`invalid digits: ${_}`);
      digits = d;
    }
    return shape;
  };
  return () => new Path(digits);
}

// node_modules/d3-shape/src/array.js
var slice = Array.prototype.slice;
function array_default(x2) {
  return typeof x2 === "object" && "length" in x2 ? x2 : Array.from(x2);
}

// node_modules/d3-shape/src/curve/linear.js
function Linear(context) {
  this._context = context;
}
Linear.prototype = {
  areaStart: function() {
    this._line = 0;
  },
  areaEnd: function() {
    this._line = NaN;
  },
  lineStart: function() {
    this._point = 0;
  },
  lineEnd: function() {
    if (this._line || this._line !== 0 && this._point === 1)
      this._context.closePath();
    this._line = 1 - this._line;
  },
  point: function(x2, y2) {
    x2 = +x2, y2 = +y2;
    switch (this._point) {
      case 0:
        this._point = 1;
        this._line ? this._context.lineTo(x2, y2) : this._context.moveTo(x2, y2);
        break;
      case 1:
        this._point = 2;
      default:
        this._context.lineTo(x2, y2);
        break;
    }
  }
};
function linear_default(context) {
  return new Linear(context);
}

// node_modules/d3-shape/src/point.js
function x(p) {
  return p[0];
}
function y(p) {
  return p[1];
}

// node_modules/d3-shape/src/line.js
function line_default(x2, y2) {
  var defined = constant_default4(true), context = null, curve = linear_default, output = null, path2 = withPath(line);
  x2 = typeof x2 === "function" ? x2 : x2 === void 0 ? x : constant_default4(x2);
  y2 = typeof y2 === "function" ? y2 : y2 === void 0 ? y : constant_default4(y2);
  function line(data) {
    var i, n = (data = array_default(data)).length, d, defined0 = false, buffer;
    if (context == null)
      output = curve(buffer = path2());
    for (i = 0; i <= n; ++i) {
      if (!(i < n && defined(d = data[i], i, data)) === defined0) {
        if (defined0 = !defined0)
          output.lineStart();
        else
          output.lineEnd();
      }
      if (defined0)
        output.point(+x2(d, i, data), +y2(d, i, data));
    }
    if (buffer)
      return output = null, buffer + "" || null;
  }
  line.x = function(_) {
    return arguments.length ? (x2 = typeof _ === "function" ? _ : constant_default4(+_), line) : x2;
  };
  line.y = function(_) {
    return arguments.length ? (y2 = typeof _ === "function" ? _ : constant_default4(+_), line) : y2;
  };
  line.defined = function(_) {
    return arguments.length ? (defined = typeof _ === "function" ? _ : constant_default4(!!_), line) : defined;
  };
  line.curve = function(_) {
    return arguments.length ? (curve = _, context != null && (output = curve(context)), line) : curve;
  };
  line.context = function(_) {
    return arguments.length ? (_ == null ? context = output = null : output = curve(context = _), line) : context;
  };
  return line;
}

// node_modules/d3-shape/src/symbol/circle.js
var circle_default = {
  draw(context, size) {
    const r = sqrt(size / pi2);
    context.moveTo(r, 0);
    context.arc(0, 0, r, 0, tau2);
  }
};

// node_modules/d3-shape/src/symbol/cross.js
var cross_default = {
  draw(context, size) {
    const r = sqrt(size / 5) / 2;
    context.moveTo(-3 * r, -r);
    context.lineTo(-r, -r);
    context.lineTo(-r, -3 * r);
    context.lineTo(r, -3 * r);
    context.lineTo(r, -r);
    context.lineTo(3 * r, -r);
    context.lineTo(3 * r, r);
    context.lineTo(r, r);
    context.lineTo(r, 3 * r);
    context.lineTo(-r, 3 * r);
    context.lineTo(-r, r);
    context.lineTo(-3 * r, r);
    context.closePath();
  }
};

// node_modules/d3-shape/src/symbol/diamond.js
var tan30 = sqrt(1 / 3);
var tan30_2 = tan30 * 2;
var diamond_default = {
  draw(context, size) {
    const y2 = sqrt(size / tan30_2);
    const x2 = y2 * tan30;
    context.moveTo(0, -y2);
    context.lineTo(x2, 0);
    context.lineTo(0, y2);
    context.lineTo(-x2, 0);
    context.closePath();
  }
};

// node_modules/d3-shape/src/symbol/triangle.js
var sqrt3 = sqrt(3);
var triangle_default = {
  draw(context, size) {
    const y2 = -sqrt(size / (sqrt3 * 3));
    context.moveTo(0, y2 * 2);
    context.lineTo(-sqrt3 * y2, -y2);
    context.lineTo(sqrt3 * y2, -y2);
    context.closePath();
  }
};

// node_modules/d3-shape/src/symbol.js
function Symbol2(type2, size) {
  let context = null, path2 = withPath(symbol);
  type2 = typeof type2 === "function" ? type2 : constant_default4(type2 || circle_default);
  size = typeof size === "function" ? size : constant_default4(size === void 0 ? 64 : +size);
  function symbol() {
    let buffer;
    if (!context)
      context = buffer = path2();
    type2.apply(this, arguments).draw(context, +size.apply(this, arguments));
    if (buffer)
      return context = null, buffer + "" || null;
  }
  symbol.type = function(_) {
    return arguments.length ? (type2 = typeof _ === "function" ? _ : constant_default4(_), symbol) : type2;
  };
  symbol.size = function(_) {
    return arguments.length ? (size = typeof _ === "function" ? _ : constant_default4(+_), symbol) : size;
  };
  symbol.context = function(_) {
    return arguments.length ? (context = _ == null ? null : _, symbol) : context;
  };
  return symbol;
}

// node_modules/d3-shape/src/noop.js
function noop_default() {
}

// node_modules/d3-shape/src/curve/cardinal.js
function point(that, x2, y2) {
  that._context.bezierCurveTo(
    that._x1 + that._k * (that._x2 - that._x0),
    that._y1 + that._k * (that._y2 - that._y0),
    that._x2 + that._k * (that._x1 - x2),
    that._y2 + that._k * (that._y1 - y2),
    that._x2,
    that._y2
  );
}
function Cardinal(context, tension) {
  this._context = context;
  this._k = (1 - tension) / 6;
}
Cardinal.prototype = {
  areaStart: function() {
    this._line = 0;
  },
  areaEnd: function() {
    this._line = NaN;
  },
  lineStart: function() {
    this._x0 = this._x1 = this._x2 = this._y0 = this._y1 = this._y2 = NaN;
    this._point = 0;
  },
  lineEnd: function() {
    switch (this._point) {
      case 2:
        this._context.lineTo(this._x2, this._y2);
        break;
      case 3:
        point(this, this._x1, this._y1);
        break;
    }
    if (this._line || this._line !== 0 && this._point === 1)
      this._context.closePath();
    this._line = 1 - this._line;
  },
  point: function(x2, y2) {
    x2 = +x2, y2 = +y2;
    switch (this._point) {
      case 0:
        this._point = 1;
        this._line ? this._context.lineTo(x2, y2) : this._context.moveTo(x2, y2);
        break;
      case 1:
        this._point = 2;
        this._x1 = x2, this._y1 = y2;
        break;
      case 2:
        this._point = 3;
      default:
        point(this, x2, y2);
        break;
    }
    this._x0 = this._x1, this._x1 = this._x2, this._x2 = x2;
    this._y0 = this._y1, this._y1 = this._y2, this._y2 = y2;
  }
};
var cardinal_default = function custom(tension) {
  function cardinal(context) {
    return new Cardinal(context, tension);
  }
  cardinal.tension = function(tension2) {
    return custom(+tension2);
  };
  return cardinal;
}(0);

// node_modules/d3-shape/src/curve/cardinalClosed.js
function CardinalClosed(context, tension) {
  this._context = context;
  this._k = (1 - tension) / 6;
}
CardinalClosed.prototype = {
  areaStart: noop_default,
  areaEnd: noop_default,
  lineStart: function() {
    this._x0 = this._x1 = this._x2 = this._x3 = this._x4 = this._x5 = this._y0 = this._y1 = this._y2 = this._y3 = this._y4 = this._y5 = NaN;
    this._point = 0;
  },
  lineEnd: function() {
    switch (this._point) {
      case 1: {
        this._context.moveTo(this._x3, this._y3);
        this._context.closePath();
        break;
      }
      case 2: {
        this._context.lineTo(this._x3, this._y3);
        this._context.closePath();
        break;
      }
      case 3: {
        this.point(this._x3, this._y3);
        this.point(this._x4, this._y4);
        this.point(this._x5, this._y5);
        break;
      }
    }
  },
  point: function(x2, y2) {
    x2 = +x2, y2 = +y2;
    switch (this._point) {
      case 0:
        this._point = 1;
        this._x3 = x2, this._y3 = y2;
        break;
      case 1:
        this._point = 2;
        this._context.moveTo(this._x4 = x2, this._y4 = y2);
        break;
      case 2:
        this._point = 3;
        this._x5 = x2, this._y5 = y2;
        break;
      default:
        point(this, x2, y2);
        break;
    }
    this._x0 = this._x1, this._x1 = this._x2, this._x2 = x2;
    this._y0 = this._y1, this._y1 = this._y2, this._y2 = y2;
  }
};
var cardinalClosed_default = function custom2(tension) {
  function cardinal(context) {
    return new CardinalClosed(context, tension);
  }
  cardinal.tension = function(tension2) {
    return custom2(+tension2);
  };
  return cardinal;
}(0);

// node_modules/d3-shape/src/curve/catmullRom.js
function point2(that, x2, y2) {
  var x1 = that._x1, y1 = that._y1, x22 = that._x2, y22 = that._y2;
  if (that._l01_a > epsilon3) {
    var a = 2 * that._l01_2a + 3 * that._l01_a * that._l12_a + that._l12_2a, n = 3 * that._l01_a * (that._l01_a + that._l12_a);
    x1 = (x1 * a - that._x0 * that._l12_2a + that._x2 * that._l01_2a) / n;
    y1 = (y1 * a - that._y0 * that._l12_2a + that._y2 * that._l01_2a) / n;
  }
  if (that._l23_a > epsilon3) {
    var b = 2 * that._l23_2a + 3 * that._l23_a * that._l12_a + that._l12_2a, m = 3 * that._l23_a * (that._l23_a + that._l12_a);
    x22 = (x22 * b + that._x1 * that._l23_2a - x2 * that._l12_2a) / m;
    y22 = (y22 * b + that._y1 * that._l23_2a - y2 * that._l12_2a) / m;
  }
  that._context.bezierCurveTo(x1, y1, x22, y22, that._x2, that._y2);
}
function CatmullRom(context, alpha) {
  this._context = context;
  this._alpha = alpha;
}
CatmullRom.prototype = {
  areaStart: function() {
    this._line = 0;
  },
  areaEnd: function() {
    this._line = NaN;
  },
  lineStart: function() {
    this._x0 = this._x1 = this._x2 = this._y0 = this._y1 = this._y2 = NaN;
    this._l01_a = this._l12_a = this._l23_a = this._l01_2a = this._l12_2a = this._l23_2a = this._point = 0;
  },
  lineEnd: function() {
    switch (this._point) {
      case 2:
        this._context.lineTo(this._x2, this._y2);
        break;
      case 3:
        this.point(this._x2, this._y2);
        break;
    }
    if (this._line || this._line !== 0 && this._point === 1)
      this._context.closePath();
    this._line = 1 - this._line;
  },
  point: function(x2, y2) {
    x2 = +x2, y2 = +y2;
    if (this._point) {
      var x23 = this._x2 - x2, y23 = this._y2 - y2;
      this._l23_a = Math.sqrt(this._l23_2a = Math.pow(x23 * x23 + y23 * y23, this._alpha));
    }
    switch (this._point) {
      case 0:
        this._point = 1;
        this._line ? this._context.lineTo(x2, y2) : this._context.moveTo(x2, y2);
        break;
      case 1:
        this._point = 2;
        break;
      case 2:
        this._point = 3;
      default:
        point2(this, x2, y2);
        break;
    }
    this._l01_a = this._l12_a, this._l12_a = this._l23_a;
    this._l01_2a = this._l12_2a, this._l12_2a = this._l23_2a;
    this._x0 = this._x1, this._x1 = this._x2, this._x2 = x2;
    this._y0 = this._y1, this._y1 = this._y2, this._y2 = y2;
  }
};
var catmullRom_default = function custom3(alpha) {
  function catmullRom(context) {
    return alpha ? new CatmullRom(context, alpha) : new Cardinal(context, 0);
  }
  catmullRom.alpha = function(alpha2) {
    return custom3(+alpha2);
  };
  return catmullRom;
}(0.5);

// node_modules/d3-shape/src/curve/catmullRomClosed.js
function CatmullRomClosed(context, alpha) {
  this._context = context;
  this._alpha = alpha;
}
CatmullRomClosed.prototype = {
  areaStart: noop_default,
  areaEnd: noop_default,
  lineStart: function() {
    this._x0 = this._x1 = this._x2 = this._x3 = this._x4 = this._x5 = this._y0 = this._y1 = this._y2 = this._y3 = this._y4 = this._y5 = NaN;
    this._l01_a = this._l12_a = this._l23_a = this._l01_2a = this._l12_2a = this._l23_2a = this._point = 0;
  },
  lineEnd: function() {
    switch (this._point) {
      case 1: {
        this._context.moveTo(this._x3, this._y3);
        this._context.closePath();
        break;
      }
      case 2: {
        this._context.lineTo(this._x3, this._y3);
        this._context.closePath();
        break;
      }
      case 3: {
        this.point(this._x3, this._y3);
        this.point(this._x4, this._y4);
        this.point(this._x5, this._y5);
        break;
      }
    }
  },
  point: function(x2, y2) {
    x2 = +x2, y2 = +y2;
    if (this._point) {
      var x23 = this._x2 - x2, y23 = this._y2 - y2;
      this._l23_a = Math.sqrt(this._l23_2a = Math.pow(x23 * x23 + y23 * y23, this._alpha));
    }
    switch (this._point) {
      case 0:
        this._point = 1;
        this._x3 = x2, this._y3 = y2;
        break;
      case 1:
        this._point = 2;
        this._context.moveTo(this._x4 = x2, this._y4 = y2);
        break;
      case 2:
        this._point = 3;
        this._x5 = x2, this._y5 = y2;
        break;
      default:
        point2(this, x2, y2);
        break;
    }
    this._l01_a = this._l12_a, this._l12_a = this._l23_a;
    this._l01_2a = this._l12_2a, this._l12_2a = this._l23_2a;
    this._x0 = this._x1, this._x1 = this._x2, this._x2 = x2;
    this._y0 = this._y1, this._y1 = this._y2, this._y2 = y2;
  }
};
var catmullRomClosed_default = function custom4(alpha) {
  function catmullRom(context) {
    return alpha ? new CatmullRomClosed(context, alpha) : new CardinalClosed(context, 0);
  }
  catmullRom.alpha = function(alpha2) {
    return custom4(+alpha2);
  };
  return catmullRom;
}(0.5);

// node_modules/d3-zoom/src/constant.js
var constant_default5 = (x2) => () => x2;

// node_modules/d3-zoom/src/event.js
function ZoomEvent(type2, {
  sourceEvent,
  target,
  transform: transform2,
  dispatch: dispatch2
}) {
  Object.defineProperties(this, {
    type: { value: type2, enumerable: true, configurable: true },
    sourceEvent: { value: sourceEvent, enumerable: true, configurable: true },
    target: { value: target, enumerable: true, configurable: true },
    transform: { value: transform2, enumerable: true, configurable: true },
    _: { value: dispatch2 }
  });
}

// node_modules/d3-zoom/src/transform.js
function Transform(k, x2, y2) {
  this.k = k;
  this.x = x2;
  this.y = y2;
}
Transform.prototype = {
  constructor: Transform,
  scale: function(k) {
    return k === 1 ? this : new Transform(this.k * k, this.x, this.y);
  },
  translate: function(x2, y2) {
    return x2 === 0 & y2 === 0 ? this : new Transform(this.k, this.x + this.k * x2, this.y + this.k * y2);
  },
  apply: function(point3) {
    return [point3[0] * this.k + this.x, point3[1] * this.k + this.y];
  },
  applyX: function(x2) {
    return x2 * this.k + this.x;
  },
  applyY: function(y2) {
    return y2 * this.k + this.y;
  },
  invert: function(location) {
    return [(location[0] - this.x) / this.k, (location[1] - this.y) / this.k];
  },
  invertX: function(x2) {
    return (x2 - this.x) / this.k;
  },
  invertY: function(y2) {
    return (y2 - this.y) / this.k;
  },
  rescaleX: function(x2) {
    return x2.copy().domain(x2.range().map(this.invertX, this).map(x2.invert, x2));
  },
  rescaleY: function(y2) {
    return y2.copy().domain(y2.range().map(this.invertY, this).map(y2.invert, y2));
  },
  toString: function() {
    return "translate(" + this.x + "," + this.y + ") scale(" + this.k + ")";
  }
};
var identity3 = new Transform(1, 0, 0);
transform.prototype = Transform.prototype;
function transform(node) {
  while (!node.__zoom)
    if (!(node = node.parentNode))
      return identity3;
  return node.__zoom;
}

// node_modules/d3-zoom/src/noevent.js
function nopropagation2(event) {
  event.stopImmediatePropagation();
}
function noevent_default3(event) {
  event.preventDefault();
  event.stopImmediatePropagation();
}

// node_modules/d3-zoom/src/zoom.js
function defaultFilter(event) {
  return (!event.ctrlKey || event.type === "wheel") && !event.button;
}
function defaultExtent() {
  var e = this;
  if (e instanceof SVGElement) {
    e = e.ownerSVGElement || e;
    if (e.hasAttribute("viewBox")) {
      e = e.viewBox.baseVal;
      return [[e.x, e.y], [e.x + e.width, e.y + e.height]];
    }
    return [[0, 0], [e.width.baseVal.value, e.height.baseVal.value]];
  }
  return [[0, 0], [e.clientWidth, e.clientHeight]];
}
function defaultTransform() {
  return this.__zoom || identity3;
}
function defaultWheelDelta(event) {
  return -event.deltaY * (event.deltaMode === 1 ? 0.05 : event.deltaMode ? 1 : 2e-3) * (event.ctrlKey ? 10 : 1);
}
function defaultTouchable() {
  return navigator.maxTouchPoints || "ontouchstart" in this;
}
function defaultConstrain(transform2, extent, translateExtent) {
  var dx0 = transform2.invertX(extent[0][0]) - translateExtent[0][0], dx1 = transform2.invertX(extent[1][0]) - translateExtent[1][0], dy0 = transform2.invertY(extent[0][1]) - translateExtent[0][1], dy1 = transform2.invertY(extent[1][1]) - translateExtent[1][1];
  return transform2.translate(
    dx1 > dx0 ? (dx0 + dx1) / 2 : Math.min(0, dx0) || Math.max(0, dx1),
    dy1 > dy0 ? (dy0 + dy1) / 2 : Math.min(0, dy0) || Math.max(0, dy1)
  );
}
function zoom_default2() {
  var filter2 = defaultFilter, extent = defaultExtent, constrain = defaultConstrain, wheelDelta = defaultWheelDelta, touchable = defaultTouchable, scaleExtent = [0, Infinity], translateExtent = [[-Infinity, -Infinity], [Infinity, Infinity]], duration = 250, interpolate = zoom_default, listeners = dispatch_default("start", "zoom", "end"), touchstarting, touchfirst, touchending, touchDelay = 500, wheelDelay = 150, clickDistance2 = 0, tapDistance = 10;
  function zoom(selection2) {
    selection2.property("__zoom", defaultTransform).on("wheel.zoom", wheeled, { passive: false }).on("mousedown.zoom", mousedowned).on("dblclick.zoom", dblclicked).filter(touchable).on("touchstart.zoom", touchstarted).on("touchmove.zoom", touchmoved).on("touchend.zoom touchcancel.zoom", touchended).style("-webkit-tap-highlight-color", "rgba(0,0,0,0)");
  }
  zoom.transform = function(collection, transform2, point3, event) {
    var selection2 = collection.selection ? collection.selection() : collection;
    selection2.property("__zoom", defaultTransform);
    if (collection !== selection2) {
      schedule(collection, transform2, point3, event);
    } else {
      selection2.interrupt().each(function() {
        gesture(this, arguments).event(event).start().zoom(null, typeof transform2 === "function" ? transform2.apply(this, arguments) : transform2).end();
      });
    }
  };
  zoom.scaleBy = function(selection2, k, p, event) {
    zoom.scaleTo(selection2, function() {
      var k0 = this.__zoom.k, k1 = typeof k === "function" ? k.apply(this, arguments) : k;
      return k0 * k1;
    }, p, event);
  };
  zoom.scaleTo = function(selection2, k, p, event) {
    zoom.transform(selection2, function() {
      var e = extent.apply(this, arguments), t0 = this.__zoom, p0 = p == null ? centroid(e) : typeof p === "function" ? p.apply(this, arguments) : p, p1 = t0.invert(p0), k1 = typeof k === "function" ? k.apply(this, arguments) : k;
      return constrain(translate(scale(t0, k1), p0, p1), e, translateExtent);
    }, p, event);
  };
  zoom.translateBy = function(selection2, x2, y2, event) {
    zoom.transform(selection2, function() {
      return constrain(this.__zoom.translate(
        typeof x2 === "function" ? x2.apply(this, arguments) : x2,
        typeof y2 === "function" ? y2.apply(this, arguments) : y2
      ), extent.apply(this, arguments), translateExtent);
    }, null, event);
  };
  zoom.translateTo = function(selection2, x2, y2, p, event) {
    zoom.transform(selection2, function() {
      var e = extent.apply(this, arguments), t = this.__zoom, p0 = p == null ? centroid(e) : typeof p === "function" ? p.apply(this, arguments) : p;
      return constrain(identity3.translate(p0[0], p0[1]).scale(t.k).translate(
        typeof x2 === "function" ? -x2.apply(this, arguments) : -x2,
        typeof y2 === "function" ? -y2.apply(this, arguments) : -y2
      ), e, translateExtent);
    }, p, event);
  };
  function scale(transform2, k) {
    k = Math.max(scaleExtent[0], Math.min(scaleExtent[1], k));
    return k === transform2.k ? transform2 : new Transform(k, transform2.x, transform2.y);
  }
  function translate(transform2, p0, p1) {
    var x2 = p0[0] - p1[0] * transform2.k, y2 = p0[1] - p1[1] * transform2.k;
    return x2 === transform2.x && y2 === transform2.y ? transform2 : new Transform(transform2.k, x2, y2);
  }
  function centroid(extent2) {
    return [(+extent2[0][0] + +extent2[1][0]) / 2, (+extent2[0][1] + +extent2[1][1]) / 2];
  }
  function schedule(transition2, transform2, point3, event) {
    transition2.on("start.zoom", function() {
      gesture(this, arguments).event(event).start();
    }).on("interrupt.zoom end.zoom", function() {
      gesture(this, arguments).event(event).end();
    }).tween("zoom", function() {
      var that = this, args = arguments, g = gesture(that, args).event(event), e = extent.apply(that, args), p = point3 == null ? centroid(e) : typeof point3 === "function" ? point3.apply(that, args) : point3, w = Math.max(e[1][0] - e[0][0], e[1][1] - e[0][1]), a = that.__zoom, b = typeof transform2 === "function" ? transform2.apply(that, args) : transform2, i = interpolate(a.invert(p).concat(w / a.k), b.invert(p).concat(w / b.k));
      return function(t) {
        if (t === 1)
          t = b;
        else {
          var l = i(t), k = w / l[2];
          t = new Transform(k, p[0] - l[0] * k, p[1] - l[1] * k);
        }
        g.zoom(null, t);
      };
    });
  }
  function gesture(that, args, clean) {
    return !clean && that.__zooming || new Gesture(that, args);
  }
  function Gesture(that, args) {
    this.that = that;
    this.args = args;
    this.active = 0;
    this.sourceEvent = null;
    this.extent = extent.apply(that, args);
    this.taps = 0;
  }
  Gesture.prototype = {
    event: function(event) {
      if (event)
        this.sourceEvent = event;
      return this;
    },
    start: function() {
      if (++this.active === 1) {
        this.that.__zooming = this;
        this.emit("start");
      }
      return this;
    },
    zoom: function(key, transform2) {
      if (this.mouse && key !== "mouse")
        this.mouse[1] = transform2.invert(this.mouse[0]);
      if (this.touch0 && key !== "touch")
        this.touch0[1] = transform2.invert(this.touch0[0]);
      if (this.touch1 && key !== "touch")
        this.touch1[1] = transform2.invert(this.touch1[0]);
      this.that.__zoom = transform2;
      this.emit("zoom");
      return this;
    },
    end: function() {
      if (--this.active === 0) {
        delete this.that.__zooming;
        this.emit("end");
      }
      return this;
    },
    emit: function(type2) {
      var d = select_default2(this.that).datum();
      listeners.call(
        type2,
        this.that,
        new ZoomEvent(type2, {
          sourceEvent: this.sourceEvent,
          target: zoom,
          type: type2,
          transform: this.that.__zoom,
          dispatch: listeners
        }),
        d
      );
    }
  };
  function wheeled(event, ...args) {
    if (!filter2.apply(this, arguments))
      return;
    var g = gesture(this, args).event(event), t = this.__zoom, k = Math.max(scaleExtent[0], Math.min(scaleExtent[1], t.k * Math.pow(2, wheelDelta.apply(this, arguments)))), p = pointer_default(event);
    if (g.wheel) {
      if (g.mouse[0][0] !== p[0] || g.mouse[0][1] !== p[1]) {
        g.mouse[1] = t.invert(g.mouse[0] = p);
      }
      clearTimeout(g.wheel);
    } else if (t.k === k)
      return;
    else {
      g.mouse = [p, t.invert(p)];
      interrupt_default(this);
      g.start();
    }
    noevent_default3(event);
    g.wheel = setTimeout(wheelidled, wheelDelay);
    g.zoom("mouse", constrain(translate(scale(t, k), g.mouse[0], g.mouse[1]), g.extent, translateExtent));
    function wheelidled() {
      g.wheel = null;
      g.end();
    }
  }
  function mousedowned(event, ...args) {
    if (touchending || !filter2.apply(this, arguments))
      return;
    var currentTarget = event.currentTarget, g = gesture(this, args, true).event(event), v = select_default2(event.view).on("mousemove.zoom", mousemoved, true).on("mouseup.zoom", mouseupped, true), p = pointer_default(event, currentTarget), x0 = event.clientX, y0 = event.clientY;
    nodrag_default(event.view);
    nopropagation2(event);
    g.mouse = [p, this.__zoom.invert(p)];
    interrupt_default(this);
    g.start();
    function mousemoved(event2) {
      noevent_default3(event2);
      if (!g.moved) {
        var dx = event2.clientX - x0, dy = event2.clientY - y0;
        g.moved = dx * dx + dy * dy > clickDistance2;
      }
      g.event(event2).zoom("mouse", constrain(translate(g.that.__zoom, g.mouse[0] = pointer_default(event2, currentTarget), g.mouse[1]), g.extent, translateExtent));
    }
    function mouseupped(event2) {
      v.on("mousemove.zoom mouseup.zoom", null);
      yesdrag(event2.view, g.moved);
      noevent_default3(event2);
      g.event(event2).end();
    }
  }
  function dblclicked(event, ...args) {
    if (!filter2.apply(this, arguments))
      return;
    var t0 = this.__zoom, p0 = pointer_default(event.changedTouches ? event.changedTouches[0] : event, this), p1 = t0.invert(p0), k1 = t0.k * (event.shiftKey ? 0.5 : 2), t1 = constrain(translate(scale(t0, k1), p0, p1), extent.apply(this, args), translateExtent);
    noevent_default3(event);
    if (duration > 0)
      select_default2(this).transition().duration(duration).call(schedule, t1, p0, event);
    else
      select_default2(this).call(zoom.transform, t1, p0, event);
  }
  function touchstarted(event, ...args) {
    if (!filter2.apply(this, arguments))
      return;
    var touches = event.touches, n = touches.length, g = gesture(this, args, event.changedTouches.length === n).event(event), started, i, t, p;
    nopropagation2(event);
    for (i = 0; i < n; ++i) {
      t = touches[i], p = pointer_default(t, this);
      p = [p, this.__zoom.invert(p), t.identifier];
      if (!g.touch0)
        g.touch0 = p, started = true, g.taps = 1 + !!touchstarting;
      else if (!g.touch1 && g.touch0[2] !== p[2])
        g.touch1 = p, g.taps = 0;
    }
    if (touchstarting)
      touchstarting = clearTimeout(touchstarting);
    if (started) {
      if (g.taps < 2)
        touchfirst = p[0], touchstarting = setTimeout(function() {
          touchstarting = null;
        }, touchDelay);
      interrupt_default(this);
      g.start();
    }
  }
  function touchmoved(event, ...args) {
    if (!this.__zooming)
      return;
    var g = gesture(this, args).event(event), touches = event.changedTouches, n = touches.length, i, t, p, l;
    noevent_default3(event);
    for (i = 0; i < n; ++i) {
      t = touches[i], p = pointer_default(t, this);
      if (g.touch0 && g.touch0[2] === t.identifier)
        g.touch0[0] = p;
      else if (g.touch1 && g.touch1[2] === t.identifier)
        g.touch1[0] = p;
    }
    t = g.that.__zoom;
    if (g.touch1) {
      var p0 = g.touch0[0], l0 = g.touch0[1], p1 = g.touch1[0], l1 = g.touch1[1], dp = (dp = p1[0] - p0[0]) * dp + (dp = p1[1] - p0[1]) * dp, dl = (dl = l1[0] - l0[0]) * dl + (dl = l1[1] - l0[1]) * dl;
      t = scale(t, Math.sqrt(dp / dl));
      p = [(p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2];
      l = [(l0[0] + l1[0]) / 2, (l0[1] + l1[1]) / 2];
    } else if (g.touch0)
      p = g.touch0[0], l = g.touch0[1];
    else
      return;
    g.zoom("touch", constrain(translate(t, p, l), g.extent, translateExtent));
  }
  function touchended(event, ...args) {
    if (!this.__zooming)
      return;
    var g = gesture(this, args).event(event), touches = event.changedTouches, n = touches.length, i, t;
    nopropagation2(event);
    if (touchending)
      clearTimeout(touchending);
    touchending = setTimeout(function() {
      touchending = null;
    }, touchDelay);
    for (i = 0; i < n; ++i) {
      t = touches[i];
      if (g.touch0 && g.touch0[2] === t.identifier)
        delete g.touch0;
      else if (g.touch1 && g.touch1[2] === t.identifier)
        delete g.touch1;
    }
    if (g.touch1 && !g.touch0)
      g.touch0 = g.touch1, delete g.touch1;
    if (g.touch0)
      g.touch0[1] = this.__zoom.invert(g.touch0[0]);
    else {
      g.end();
      if (g.taps === 2) {
        t = pointer_default(t, this);
        if (Math.hypot(touchfirst[0] - t[0], touchfirst[1] - t[1]) < tapDistance) {
          var p = select_default2(this).on("dblclick.zoom");
          if (p)
            p.apply(this, arguments);
        }
      }
    }
  }
  zoom.wheelDelta = function(_) {
    return arguments.length ? (wheelDelta = typeof _ === "function" ? _ : constant_default5(+_), zoom) : wheelDelta;
  };
  zoom.filter = function(_) {
    return arguments.length ? (filter2 = typeof _ === "function" ? _ : constant_default5(!!_), zoom) : filter2;
  };
  zoom.touchable = function(_) {
    return arguments.length ? (touchable = typeof _ === "function" ? _ : constant_default5(!!_), zoom) : touchable;
  };
  zoom.extent = function(_) {
    return arguments.length ? (extent = typeof _ === "function" ? _ : constant_default5([[+_[0][0], +_[0][1]], [+_[1][0], +_[1][1]]]), zoom) : extent;
  };
  zoom.scaleExtent = function(_) {
    return arguments.length ? (scaleExtent[0] = +_[0], scaleExtent[1] = +_[1], zoom) : [scaleExtent[0], scaleExtent[1]];
  };
  zoom.translateExtent = function(_) {
    return arguments.length ? (translateExtent[0][0] = +_[0][0], translateExtent[1][0] = +_[1][0], translateExtent[0][1] = +_[0][1], translateExtent[1][1] = +_[1][1], zoom) : [[translateExtent[0][0], translateExtent[0][1]], [translateExtent[1][0], translateExtent[1][1]]];
  };
  zoom.constrain = function(_) {
    return arguments.length ? (constrain = _, zoom) : constrain;
  };
  zoom.duration = function(_) {
    return arguments.length ? (duration = +_, zoom) : duration;
  };
  zoom.interpolate = function(_) {
    return arguments.length ? (interpolate = _, zoom) : interpolate;
  };
  zoom.on = function() {
    var value = listeners.on.apply(listeners, arguments);
    return value === listeners ? zoom : value;
  };
  zoom.clickDistance = function(_) {
    return arguments.length ? (clickDistance2 = (_ = +_) * _, zoom) : Math.sqrt(clickDistance2);
  };
  zoom.tapDistance = function(_) {
    return arguments.length ? (tapDistance = +_, zoom) : tapDistance;
  };
  return zoom;
}

// js/flux_network/rendering.ts
var CSS = `
.flov-wrapper{background:#0d1117;position:relative;font-family:'Inter',Arial,sans-serif;overflow:hidden;}
.flov-toolbar{display:flex;flex-wrap:wrap;gap:6px;padding:8px 12px;background:#161b22;border-bottom:1px solid #30363d;align-items:flex-end;}
.flov-btn-group{display:flex;gap:3px;align-items:flex-end;}
.flov-btn-label{font:10px/1 'Inter',Arial,sans-serif;color:#7d8590;padding:0 4px 5px 0;}
.flov-btn{background:#21262d;color:#e6edf3;border:1px solid #30363d;border-radius:6px;padding:5px 13px;font:11px 'Inter',Arial,sans-serif;cursor:pointer;transition:background 80ms,border-color 80ms;white-space:nowrap;}
.flov-btn:hover{background:#30363d;border-color:#8b949e;}
.flov-btn.active{background:#1f6feb;border-color:#1f6feb;color:#fff;}
.flov-svg{display:block;}
.flov-tt{position:absolute;pointer-events:none;z-index:1000;background:#1c2128;border-radius:8px;padding:10px 14px;font:12px/1.65 'Inter',Arial,sans-serif;color:#e6edf3;max-width:290px;box-shadow:0 6px 24px rgba(0,0,0,0.75);display:none;border:1.5px solid #30363d;}
.flov-tt b{color:#fff;}
.flov-tt .ft-flux{font-family:monospace;color:#ffa657;}
.flov-tt .ft-div{color:#30363d;}
.flov-tt .ft-comp{color:#8b949e;font-style:italic;}
.flov-tt .ft-badge{background:#21262d;border-radius:3px;padding:1px 5px;font-size:10px;color:#8b949e;letter-spacing:.04em;}
.flov-search-wrap{position:absolute;top:10px;right:10px;z-index:100;background:#161b22;border:1px solid #30363d;border-radius:6px;padding:5px 10px;width:230px;box-shadow:0 4px 14px rgba(0,0,0,0.55);}
.flov-search-input{background:transparent;border:none;outline:none;color:#e6edf3;width:100%;font:11px 'Inter',Arial,sans-serif;}
.flov-search-input::placeholder{color:#7d8590;}
.flov-search-results{display:none;max-height:184px;overflow-y:auto;margin-top:4px;}
.flov-search-result{padding:3px 8px;cursor:pointer;font-size:11px;border-top:1px solid #21262d;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;color:#e6edf3;}
.flov-search-result:hover{background:#21262d;}
.flov-station-panel{position:absolute;top:52px;left:10px;z-index:95;width:min(720px,calc(100% - 20px));max-height:min(520px,calc(100% - 70px));overflow-y:auto;background:#161b22;border:1px solid #30363d;border-radius:8px;box-shadow:0 10px 28px rgba(0,0,0,0.62);color:#e6edf3;font:12px/1.45 'Inter',Arial,sans-serif;display:none;}
.flov-station-head{position:sticky;top:0;background:#161b22;border-bottom:1px solid #30363d;padding:10px 12px;z-index:1;}
.flov-station-title{display:flex;align-items:center;gap:8px;font-weight:700;color:#fff;}
.flov-station-badge{background:#21262d;border-radius:3px;padding:1px 5px;font-size:10px;color:#8b949e;letter-spacing:.04em;}
.flov-station-close{margin-left:auto;background:transparent;border:0;color:#8b949e;font:18px/1 Arial,sans-serif;cursor:pointer;padding:0 2px;}
.flov-station-close:hover{color:#e6edf3;}
.flov-station-sub{margin-top:3px;color:#8b949e;font-size:11px;}
.flov-station-body{padding:7px 8px 10px;}
.flov-station-group{border-top:1px solid #21262d;}
.flov-station-group:first-child{border-top:0;}
.flov-station-toggle{display:flex;align-items:center;gap:7px;width:100%;padding:8px 5px;background:transparent;border:0;color:#e6edf3;text-align:left;font:12px/1.2 'Inter',Arial,sans-serif;cursor:pointer;}
.flov-station-toggle:hover{background:#21262d;}
.flov-station-caret{width:12px;color:#8b949e;font-size:11px;}
.flov-station-swatch{width:9px;height:9px;border-radius:2px;flex:0 0 auto;}
.flov-station-count{margin-left:auto;color:#7d8590;font-size:10px;}
.flov-station-rxns{display:none;padding:0 5px 8px 31px;}
.flov-station-group.open .flov-station-rxns{display:block;}
.flov-station-map{width:100%;height:auto;margin:5px 0 10px;overflow:visible;}
.flov-station-branch{cursor:pointer;}
.flov-station-branch:hover .flov-station-map-line{stroke-width:2.4;}
.flov-station-branch.search-hit .flov-station-map-line{stroke-width:3.1;}
.flov-station-branch.search-hidden{display:none;}
.flov-station-map.search-hidden{display:none;}
.flov-station-group.search-hidden{display:none;}
.flov-station-map-line{stroke-linecap:round;stroke-linejoin:round;}
.flov-station-map-label{font:12px ui-monospace,SFMono-Regular,Menlo,monospace;fill:#8b949e;}
.flov-station-map-flux{font:10px ui-monospace,SFMono-Regular,Menlo,monospace;fill:#7d8590;}
.flov-station-map-met{font:9px 'Inter',Arial,sans-serif;fill:#8b949e;}
.flov-legend{position:absolute;top:52px;right:8px;z-index:90;width:230px;overflow-y:auto;background:#161b22;border:1px solid #30363d;border-radius:8px;box-shadow:0 4px 14px rgba(0,0,0,0.55);color:#e6edf3;font:11px/1.45 'Inter',Arial,sans-serif;}
.flov-legend-head{display:flex;align-items:center;gap:6px;padding:7px 10px;border-bottom:1px solid #30363d;cursor:pointer;user-select:none;background:#1c2128;border-top-left-radius:8px;border-top-right-radius:8px;}
.flov-legend-title{font-weight:700;color:#fff;flex:1;font-size:12px;letter-spacing:.02em;}
.flov-legend-caret{color:#8b949e;font-size:11px;}
.flov-legend-body{padding:4px 10px 10px;}
.flov-legend.collapsed .flov-legend-body{display:none;}
.flov-legend.collapsed{border-bottom-left-radius:8px;border-bottom-right-radius:8px;}
.flov-legend.collapsed .flov-legend-head{border-bottom:0;}
.flov-legend-section{margin:8px 0 4px;}
.flov-legend-section + .flov-legend-section{border-top:1px solid #21262d;padding-top:6px;}
.flov-legend-section-title{font-size:9px;letter-spacing:.08em;text-transform:uppercase;color:#8b949e;margin-bottom:4px;font-weight:600;}
.flov-legend-row{display:flex;align-items:center;gap:8px;margin:3px 0;}
.flov-legend-swatch{flex:0 0 30px;display:flex;justify-content:center;align-items:center;}
.flov-legend-text{flex:1;min-width:0;}
.flov-legend-label{color:#e6edf3;font-size:11px;line-height:1.25;}
.flov-legend-sub{color:#7d8590;font-size:9.5px;line-height:1.2;margin-top:1px;}
.flov-legend-fluxbar{margin:2px 0 0;}
.flov-legend-note{color:#7d8590;font-size:9.5px;margin-top:2px;font-style:italic;}
.flov-wrapper.light-mode .flov-legend{background:#ffffff;border-color:#d0d7de;color:#24292f;box-shadow:0 8px 24px rgba(31,35,40,0.12);}
.flov-wrapper.light-mode .flov-legend-head{background:#f6f8fa;border-bottom-color:#d0d7de;}
.flov-wrapper.light-mode .flov-legend-title{color:#24292f;}
.flov-wrapper.light-mode .flov-legend-caret{color:#57606a;}
.flov-wrapper.light-mode .flov-legend-section + .flov-legend-section{border-top-color:#d8dee4;}
.flov-wrapper.light-mode .flov-legend-section-title,.flov-wrapper.light-mode .flov-legend-sub,.flov-wrapper.light-mode .flov-legend-note{color:#57606a;}
.flov-wrapper.light-mode .flov-legend-label{color:#24292f;}
@keyframes flov-flow-forward{from{stroke-dashoffset:0}to{stroke-dashoffset:-28}}
@keyframes flov-flow-reverse{from{stroke-dashoffset:0}to{stroke-dashoffset:28}}
.flov-flow-arrow{stroke-dasharray:9 5;}
.flov-flow-forward{animation:flov-flow-forward 1s linear infinite;}
.flov-flow-reverse{animation:flov-flow-reverse 1s linear infinite;}
.flov-wrapper.light-mode{background:#f6f8fa;color:#24292f;}
.flov-wrapper.light-mode .flov-toolbar{background:#ffffff;border-bottom-color:#d0d7de;}
.flov-wrapper.light-mode .flov-btn-group{color:#24292f;}
.flov-wrapper.light-mode .flov-btn-label{color:#57606a;}
.flov-wrapper.light-mode .flov-btn{background:#f6f8fa;color:#24292f;border-color:#d0d7de;}
.flov-wrapper.light-mode .flov-btn:hover{background:#eef2f7;border-color:#8c959f;}
.flov-wrapper.light-mode .flov-btn.active{background:#0969da;border-color:#0969da;color:#ffffff;}
.flov-wrapper.light-mode .flov-tt{background:#ffffff;color:#24292f;border-color:#d0d7de;box-shadow:0 10px 28px rgba(31,35,40,0.16);}
.flov-wrapper.light-mode .flov-tt b{color:#24292f;}
.flov-wrapper.light-mode .flov-tt .ft-flux{color:#b35900;}
.flov-wrapper.light-mode .flov-tt .ft-div,.flov-wrapper.light-mode .flov-tt .ft-comp,.flov-wrapper.light-mode .flov-tt .ft-badge{color:#57606a;}
.flov-wrapper.light-mode .flov-tt .ft-badge{background:#f3f4f6;}
.flov-wrapper.light-mode .flov-search-wrap{background:#ffffff;border-color:#d0d7de;box-shadow:0 8px 24px rgba(31,35,40,0.12);}
.flov-wrapper.light-mode .flov-search-input{color:#24292f;}
.flov-wrapper.light-mode .flov-search-input::placeholder{color:#57606a;}
.flov-wrapper.light-mode .flov-search-result{color:#24292f;border-top-color:#d8dee4;}
.flov-wrapper.light-mode .flov-search-result:hover{background:#eef2f7;}
.flov-wrapper.light-mode .flov-station-panel{background:#ffffff;color:#24292f;border-color:#d0d7de;box-shadow:0 10px 28px rgba(31,35,40,0.16);}
.flov-wrapper.light-mode .flov-station-head{background:#ffffff;border-bottom-color:#d0d7de;}
.flov-wrapper.light-mode .flov-station-title{color:#24292f;}
.flov-wrapper.light-mode .flov-station-badge{background:#f3f4f6;color:#57606a;}
.flov-wrapper.light-mode .flov-station-close{color:#57606a;}
.flov-wrapper.light-mode .flov-station-close:hover{color:#24292f;}
.flov-wrapper.light-mode .flov-station-sub,.flov-wrapper.light-mode .flov-station-count{color:#57606a;}
.flov-wrapper.light-mode .flov-station-group{border-top-color:#d8dee4;}
.flov-wrapper.light-mode .flov-station-toggle{color:#24292f;}
.flov-wrapper.light-mode .flov-station-toggle:hover{background:#eef2f7;}
.flov-wrapper.light-mode .flov-station-map-label,.flov-wrapper.light-mode .flov-station-map-flux,.flov-wrapper.light-mode .flov-station-map-met{fill:#57606a;}
`;
var _cssInjected = false;
function injectCSS() {
  if (_cssInjected)
    return;
  const s = document.createElement("style");
  s.textContent = CSS;
  document.head.appendChild(s);
  _cssInjected = true;
}
var symbolHexagon = {
  draw(ctx, size) {
    const r = Math.sqrt(size / (3 * Math.sqrt(3) / 2));
    for (let i = 0; i < 6; i++) {
      const a = i * Math.PI / 3;
      i === 0 ? ctx.moveTo(r * Math.cos(a), r * Math.sin(a)) : ctx.lineTo(r * Math.cos(a), r * Math.sin(a));
    }
    ctx.closePath();
  }
};
function symArea(r) {
  return Math.max(Math.PI * (r / 2) * (r / 2), 6);
}
function symPath(symbol, r) {
  const area = symArea(r);
  switch (symbol) {
    case "triangle-up":
      return Symbol2(triangle_default, area)(null) ?? "";
    case "triangle-down":
      return Symbol2(triangle_default, area)(null) ?? "";
    case "hexagon":
      return Symbol2(symbolHexagon, area)(null) ?? "";
    case "cross":
      return Symbol2(cross_default, area * 0.7)(null) ?? "";
    default:
      return Symbol2(circle_default, area)(null) ?? "";
  }
}
function metPath(m) {
  if (m.display_kind === "stop") {
    return Symbol2(circle_default, Math.PI * 20)(null) ?? "";
  }
  return Symbol2(diamond_default, Math.PI * 16)(null) ?? "";
}
function hullPath(verts, xS, yS) {
  const gen = line_default().x((d) => xS(d[0])).y((d) => yS(d[1])).curve(catmullRomClosed_default.alpha(0.5));
  return gen(verts) ?? "";
}
var LANE_W_PX = 3.5;
var MAX_RENDER_LANE_SLOT = 6;
function renderLaneSlot(slot) {
  if (!Number.isFinite(slot))
    return 0;
  return MAX_RENDER_LANE_SLOT * Math.tanh(slot / MAX_RENDER_LANE_SLOT);
}
function railOffsetPath(pts, seg_slots, xS, yS, k_zoom) {
  const n = pts.length;
  if (n < 2)
    return "";
  const sx = pts.map((p) => xS(p[0]));
  const sy = pts.map((p) => yS(p[1]));
  const ux = [], uy = [];
  const nx = [], ny = [];
  for (let k = 0; k < n - 1; k++) {
    const dx = sx[k + 1] - sx[k], dy = sy[k + 1] - sy[k];
    const len = Math.hypot(dx, dy) || 1;
    ux.push(dx / len);
    uy.push(dy / len);
    nx.push(dy / len);
    ny.push(-dx / len);
  }
  const off = seg_slots.map((s) => renderLaneSlot(s) * LANE_W_PX / k_zoom);
  const MITER_CAP = 4 * LANE_W_PX / k_zoom;
  const pts_out = [];
  pts_out.push([sx[0] + off[0] * nx[0], sy[0] + off[0] * ny[0]]);
  for (let k = 1; k < n - 1; k++) {
    const Ax = sx[k] + off[k - 1] * nx[k - 1], Ay = sy[k] + off[k - 1] * ny[k - 1];
    const Bx = sx[k] + off[k] * nx[k], By = sy[k] + off[k] * ny[k];
    const cross = ux[k - 1] * uy[k] - uy[k - 1] * ux[k];
    if (Math.abs(cross) < 1e-9) {
      if (Math.abs(off[k] - off[k - 1]) < 1e-6) {
        pts_out.push([(Ax + Bx) / 2, (Ay + By) / 2]);
      } else {
        pts_out.push([Ax, Ay]);
        pts_out.push([Bx, By]);
      }
      continue;
    } else {
      const ddx = Bx - Ax, ddy = By - Ay;
      const t = (ddx * uy[k] - ddy * ux[k]) / cross;
      const mx = Ax + t * ux[k - 1], my = Ay + t * uy[k - 1];
      const dist = Math.hypot(mx - sx[k], my - sy[k]);
      if (dist > MITER_CAP) {
        pts_out.push([(Ax + Bx) / 2, (Ay + By) / 2]);
      } else {
        pts_out.push([mx, my]);
      }
    }
  }
  pts_out.push([sx[n - 1] + off[n - 2] * nx[n - 2], sy[n - 1] + off[n - 2] * ny[n - 2]]);
  let d = `M${pts_out[0][0].toFixed(1)} ${pts_out[0][1].toFixed(1)}`;
  for (let i = 1; i < pts_out.length; i++) {
    d += ` L${pts_out[i][0].toFixed(1)} ${pts_out[i][1].toFixed(1)}`;
  }
  return d;
}
var HTML_ESCAPES = {
  "&": "&amp;",
  "<": "&lt;",
  ">": "&gt;",
  '"': "&quot;",
  "'": "&#39;"
};
function escapeHTML(value) {
  return String(value ?? "").replace(/[&<>"']/g, (c) => HTML_ESCAPES[c]);
}
function rxnTooltipHTML(h) {
  const badge = h.kind_badge ? `<span class="ft-badge">${escapeHTML(h.kind_badge)}</span>` : "";
  const id2 = h.id ? ` <span style="color:#8b949e;font-size:11px">(${escapeHTML(h.id)})</span>` : "";
  const pipe = h.pipeline_tag ? ` \xB7 ${escapeHTML(h.pipeline_tag)}` : "";
  const std = h.std_str ? `<br><span style="color:#7d8590;font-size:10px">\u03C3 = ${escapeHTML(h.std_str)}</span>` : "";
  const extra = h.extra ? `<br>${escapeHTML(h.extra)}` : "";
  return `<b>${escapeHTML(h.display_name)}</b>${id2} ${badge}<br><span class="ft-div">\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500</span><br><span class="ft-comp">${escapeHTML(h.comp_label)}${pipe}</span><br>Flux: <span class="ft-flux">${escapeHTML(h.flux_str)}</span>${std}<br><span style="color:#8b949e">Sub:</span> ${escapeHTML(h.substrates)}<br><span style="color:#8b949e">Prd:</span> ${escapeHTML(h.products)}${extra}`;
}
function metTooltipHTML(m) {
  return `<b>${escapeHTML(m.name)}</b><br><span class="ft-div">\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500</span><br><span class="ft-comp">${escapeHTML(m.comp_label)}</span><br><span style="color:#8b949e">ID:</span> ${escapeHTML(m.id)}<br><span style="color:#8b949e">Consumed:</span> ${escapeHTML(m.consumers.join(", ") || "\u2014")}<br><span style="color:#8b949e">Produced:</span> ${escapeHTML(m.producers.join(", ") || "\u2014")}`;
}
var KLASS_COLORS = {
  amino_acid: "#22c55e",
  sugar: "#f59e0b",
  cofactor: "#a855f7",
  inorganic: "#06b6d4",
  other: "#ec4899"
};
var STATION_BRANCH_COLORS = [
  "#22c55e",
  "#f59e0b",
  "#38bdf8",
  "#ec4899",
  "#a855f7"
];

// js/flux_network/index.ts
function render({ model: S, el: I }) {
  injectCSS();
  const wrapper = document.createElement("div");
  wrapper.className = "flov-wrapper";
  I.appendChild(wrapper);
  let fig = null;
  let currentView = 0;
  let xS = linear2();
  let yS = linear2();
  let svgEl = null;
  let zoomB = null;
  let viewBtns = [];
  let tooltipEl = null;
  let stationPanelEl = null;
  let activeStationPanel = null;
  let currentStationSearch = "";
  let resizeTimer = null;
  let gMetSel = null;
  let gMetHitSel = null;
  let gRxnSel = null;
  let gRxnHitSel = null;
  let currentTransform = identity3;
  let rafId = null;
  let lightMode = true;
  let focusState = null;
  const arrowMarkerEndId = `flov-arrow-end-${Math.random().toString(36).slice(2)}`;
  const arrowMarkerStartId = `flov-arrow-start-${Math.random().toString(36).slice(2)}`;
  const MARGIN = { t: 10, r: 90, b: 10, l: 20 };
  const MET_HOVER_STROKE = 12;
  const RXN_HOVER_STROKE = 14;
  function metTransform(node, transform2) {
    const tx = xS(node.x).toFixed(1);
    const ty = yS(node.y).toFixed(1);
    const inv = (1 / transform2.k).toFixed(6);
    return `matrix(${inv},0,0,${inv},${tx},${ty})`;
  }
  function updateMetZoom(transform2) {
    if (!gMetSel)
      return;
    gMetSel.selectAll("path.met-shape").attr("transform", (d) => metTransform(d, transform2));
    gMetHitSel?.attr("transform", (d) => metTransform(d, transform2));
  }
  function updateRailZoom(transform2) {
    const zl = wrapper.querySelector(".zoom-layer");
    if (!zl || !fig)
      return;
    select_default2(zl).select(".g-rail").selectAll("path.rail-line").attr("d", (d) => railOffsetPath(d.path, d.seg_slots, xS, yS, transform2.k));
  }
  function updateRouteZoom(transform2) {
    const zl = wrapper.querySelector(".zoom-layer");
    if (!zl || !fig)
      return;
    select_default2(zl).select(".g-route").selectAll("path.route-line").attr("d", (d) => railOffsetPath(d.path, d.seg_slots, xS, yS, transform2.k));
  }
  function updateRxnZoom(transform2) {
    if (!gRxnSel || !gRxnHitSel)
      return;
    gRxnSel.attr("transform", (d) => rxnTransform(d, transform2));
    gRxnHitSel.attr("transform", (d) => rxnTransform(d, transform2));
  }
  function buildUI(data) {
    fig = data;
    currentView = 0;
    currentTransform = identity3;
    gMetSel = null;
    gMetHitSel = null;
    gRxnSel = null;
    gRxnHitSel = null;
    focusState = null;
    activeStationPanel = null;
    if (rafId !== null) {
      cancelAnimationFrame(rafId);
      rafId = null;
    }
    wrapper.innerHTML = "";
    tooltipEl = null;
    stationPanelEl = null;
    viewBtns = [];
    const svgW = Math.max(wrapper.clientWidth || 900, 600);
    const svgH = data.meta.height;
    xS = linear2().domain(data.meta.x_range).range([MARGIN.l, svgW - MARGIN.r]);
    yS = linear2().domain(data.meta.y_range).range([svgH - MARGIN.b, MARGIN.t]);
    wrapper.classList.toggle("light-mode", lightMode);
    buildToolbar(data, wrapper);
    const { svg, zoomLayer } = buildSVG(svgW, svgH);
    buildLegend(wrapper, data, svgH);
    buildSearch(wrapper, data);
    buildStationPanel(wrapper);
    buildTooltip(wrapper);
    drawAll(zoomLayer, data, svgW);
    zoomB = zoom_default2().scaleExtent([0.12, 40]).on("zoom", (e) => {
      currentTransform = e.transform;
      zoomLayer.attr("transform", e.transform.toString());
      if (rafId !== null)
        cancelAnimationFrame(rafId);
      rafId = requestAnimationFrame(() => {
        rafId = null;
        updateRxnZoom(currentTransform);
        updateMetZoom(currentTransform);
        updateRailZoom(currentTransform);
        updateRouteZoom(currentTransform);
      });
    });
    svg.call(zoomB);
    svg.on("dblclick.zoom", () => {
      svg.transition().duration(380).call(zoomB.transform, identity3);
    });
    svg.on("click.focus", () => clearFocus());
  }
  function buildSVG(svgW, svgH) {
    const svg = select_default2(wrapper).append("svg").attr("class", "flov-svg").attr("width", svgW).attr("height", svgH).style("background", lightMode ? "#f6f8fa" : "#0d1117");
    const defs = svg.append("defs");
    defs.append("marker").attr("id", arrowMarkerEndId).attr("viewBox", "0 0 8 8").attr("refX", 7.2).attr("refY", 4).attr("markerWidth", 3.5).attr("markerHeight", 3.5).attr("orient", "auto").append("path").attr("d", "M0,0 L8,4 L0,8 Z").attr("fill", "context-stroke");
    defs.append("marker").attr("id", arrowMarkerStartId).attr("viewBox", "0 0 8 8").attr("refX", 7.2).attr("refY", 4).attr("markerWidth", 3.5).attr("markerHeight", 3.5).attr("orient", "auto-start-reverse").append("path").attr("d", "M0,0 L8,4 L0,8 Z").attr("fill", "context-stroke");
    svgEl = svg.node();
    const zoomLayer = svg.append("g").attr("class", "zoom-layer").style("will-change", "transform");
    return { svg, zoomLayer };
  }
  function drawAll(zoomLayer, data, svgW) {
    const gComp = zoomLayer.append("g").attr("class", "g-comp");
    const gGuide = zoomLayer.append("g").attr("class", "g-guide");
    const gInner = zoomLayer.append("g").attr("class", "g-inner").style("display", "none");
    const gRail = zoomLayer.append("g").attr("class", "g-rail");
    const gRoute = zoomLayer.append("g").attr("class", "g-route");
    const gStn = zoomLayer.append("g").attr("class", "g-stn");
    gMetSel = zoomLayer.append("g").attr("class", "g-met");
    const gRxn = zoomLayer.append("g").attr("class", "g-rxn");
    const gLbl = zoomLayer.append("g").attr("class", "g-lbl").style("display", "none");
    const view0 = data.views[0];
    gComp.selectAll("path").data(data.compartments).join("path").attr("d", (c) => hullPath(c.hull_vertices, xS, yS)).attr("fill", (c) => c.fill).attr("stroke", (c) => c.color).attr("stroke-width", 1.8).attr("stroke-dasharray", "5,3").attr("vector-effect", "non-scaling-stroke").style("pointer-events", "none");
    gComp.selectAll("text").data(data.compartments).join("text").attr("x", (c) => xS(c.label_x)).attr("y", (c) => yS(c.label_y)).attr("text-anchor", "middle").attr("dominant-baseline", "auto").attr("fill", (c) => c.color).attr("font-size", 12).attr("font-weight", "bold").attr("font-family", "'Inter',Arial,sans-serif").style("pointer-events", "none").text((c) => c.label);
    gGuide.selectAll("path").data(data.guide_links).join("path").attr("d", (g) => `M${xS(g.x[0]).toFixed(1)} ${yS(g.y[0]).toFixed(1)} L${xS(g.x[1]).toFixed(1)} ${yS(g.y[1]).toFixed(1)}`).attr("stroke", "rgba(139,148,158,0.55)").attr("stroke-width", 0.9).attr("fill", "none").attr("stroke-dasharray", "4,3").attr("vector-effect", "non-scaling-stroke").style("pointer-events", "none");
    const railStyleById = /* @__PURE__ */ new Map();
    view0.rail_routes.forEach((rt) => railStyleById.set(rt.line_id, rt));
    gRail.selectAll("path").data(data.rail_routes).join("path").attr("class", "rail-line").attr("d", (d) => railOffsetPath(d.path, d.seg_slots, xS, yS, currentTransform.k)).attr("stroke", (d) => railStyleById.get(d.line_id)?.color ?? "rgba(110,118,129,0.18)").attr("stroke-width", (d) => railStyleById.get(d.line_id)?.width ?? 0.8).attr("fill", "none").attr("stroke-linecap", "round").attr("stroke-linejoin", "round").attr("vector-effect", "non-scaling-stroke").style("pointer-events", "none");
    gMetHitSel = gMetSel.selectAll("path.met-hit").data(data.met_nodes).join("path").attr("class", "met-hit").attr("d", (d) => metPath(d)).attr("transform", (d) => metTransform(d, currentTransform)).attr("fill", "transparent").attr("stroke", "transparent").attr("stroke-width", MET_HOVER_STROKE).attr("vector-effect", "non-scaling-stroke").style("pointer-events", "all").style("cursor", "pointer").on("mouseover", (ev, d) => showTT(ev, metTooltipHTML(d), d.comp_color)).on("mousemove", moveTT).on("mouseout", hideTT).on("click", (ev, d) => setMetFocus(d, ev));
    gMetSel.selectAll("path.met-shape").data(data.met_nodes).join("path").attr("class", "met-shape").attr("d", (d) => metPath(d)).attr("transform", (d) => metTransform(d, currentTransform)).attr("fill", (d) => d.display_kind === "stop" ? "#ffffff" : d.comp_color).attr("fill-opacity", (d) => d.display_kind === "stop" ? 1 : 0.85).attr("stroke", (d) => d.display_kind === "stop" ? d.comp_color : "#0d1117").attr("stroke-width", (d) => d.display_kind === "stop" ? 1.3 : 0.7).style("cursor", "pointer").on("mouseover", (ev, d) => showTT(ev, metTooltipHTML(d), d.comp_color)).on("mousemove", moveTT).on("mouseout", hideTT).on("click", (ev, d) => setMetFocus(d, ev));
    const routeData = [];
    data.stations.forEach((st) => {
      st.lines.forEach((ln, li) => {
        routeData.push({
          pair_id: st.pair_id,
          klass: ln.klass,
          rxn_ids: ln.rxn_ids,
          seg_slots: ln.seg_slots,
          path: ln.path,
          lineIndex: li
        });
      });
    });
    function lineStyle(d, view) {
      const sv = view.stations.find((s) => s.pair_id === d.pair_id);
      if (!sv || !sv.lines[d.lineIndex]) {
        return { color: "#6e7681", width: 2 };
      }
      const lv = sv.lines[d.lineIndex];
      return { color: lv.color, width: lv.width };
    }
    gRoute.selectAll("path").data(routeData).join("path").attr("class", "route-line").attr("d", (d) => railOffsetPath(d.path, d.seg_slots, xS, yS, currentTransform.k)).attr("fill", "none").attr("stroke-linejoin", "round").attr("stroke-linecap", "round").attr("vector-effect", "non-scaling-stroke").attr("stroke", (d) => lineStyle(d, view0).color).attr("stroke-width", (d) => lineStyle(d, view0).width).style("pointer-events", "none");
    const pillData = [];
    data.stations.forEach((st) => {
      pillData.push({ st, side: "a" });
      pillData.push({ st, side: "b" });
    });
    function pillRotation(p) {
      const t = p.side === "a" ? p.st.tangent_a : p.st.tangent_b;
      const ang = Math.atan2(-t[1], t[0]) * 180 / Math.PI + 90;
      return ang;
    }
    function pillTransform(p) {
      const a = p.side === "a" ? p.st.anchor_a : p.st.anchor_b;
      return `translate(${xS(a[0]).toFixed(1)},${yS(a[1]).toFixed(1)}) rotate(${pillRotation(p).toFixed(1)})`;
    }
    function pillPixelHeight(p) {
      const px0 = yS(0), px1 = yS(p.st.pill_height);
      return Math.max(7, Math.abs(px1 - px0));
    }
    function pillPixelWidth(p) {
      return Math.max(3, pillPixelHeight(p) * 0.35);
    }
    gStn.selectAll("rect.stn-pill").data(pillData).join("rect").attr("class", "stn-pill").attr("x", (p) => -pillPixelWidth(p) / 2).attr("y", (p) => -pillPixelHeight(p) / 2).attr("width", (p) => pillPixelWidth(p)).attr("height", (p) => pillPixelHeight(p)).attr("rx", (p) => pillPixelWidth(p) / 2).attr("ry", (p) => pillPixelWidth(p) / 2).attr("fill", "#ffffff").attr("stroke", (p) => p.side === "a" ? p.st.comp_a_color : p.st.comp_b_color).attr("stroke-width", 2.2).attr("vector-effect", "non-scaling-stroke").attr("transform", pillTransform).style("cursor", "pointer").on("click", (ev, p) => {
      setStationFocus(p.st, ev);
      showStationPanel(p.st);
    });
    gInner.selectAll("path").data(data.inner_edges).join("path").attr("d", (e) => `M${xS(e.x[0]).toFixed(1)} ${yS(e.y[0]).toFixed(1)} L${xS(e.x[1]).toFixed(1)} ${yS(e.y[1]).toFixed(1)}`).attr("stroke", "rgba(139,148,158,0.50)").attr("stroke-width", 0.8).attr("stroke-dasharray", "3,3").attr("fill", "none").attr("vector-effect", "non-scaling-stroke").style("pointer-events", "none");
    gRxnHitSel = gRxn.selectAll("path.rxn-hit").data(view0.rxn_nodes).join("path").attr("class", "rxn-hit").attr("d", (d) => symPath(d.symbol, d.r)).attr("transform", (d) => rxnTransform(d)).attr("fill", "transparent").attr("stroke", "transparent").attr("stroke-width", RXN_HOVER_STROKE).attr("vector-effect", "non-scaling-stroke").style("pointer-events", "all").style("cursor", "pointer").on("mouseover", (ev, d) => showTT(ev, rxnTooltipHTML(d.hover), d.border_color)).on("mousemove", moveTT).on("mouseout", hideTT).on("click", (ev, d) => setRxnFocus(d, ev));
    gRxnSel = gRxn.selectAll("path.rxn-shape").data(view0.rxn_nodes).join("path").attr("class", "rxn-shape").attr("d", (d) => symPath(d.symbol, d.r)).attr("transform", (d) => rxnTransform(d)).attr("fill", (d) => d.fill_color).attr("stroke", (d) => d.border_color).attr("stroke-width", (d) => d.border_width).attr("opacity", (d) => d.opacity).style("cursor", "pointer").on("mouseover", (ev, d) => showTT(ev, rxnTooltipHTML(d.hover), d.border_color)).on("mousemove", moveTT).on("mouseout", hideTT).on("click", (ev, d) => setRxnFocus(d, ev));
    gLbl.selectAll("text").data(view0.rxn_nodes).join("text").attr("x", (d) => xS(d.x)).attr("y", (d) => yS(d.y) - d.r - 2).attr("text-anchor", "middle").attr("fill", "#8b949e").attr("font-size", 7).attr("font-family", "'Inter',Arial,sans-serif").style("pointer-events", "none").text((d) => d.name);
    applyFocus();
  }
  function rxnTransform(d, transform2 = currentTransform) {
    const tx = xS(d.x).toFixed(1), ty = yS(d.y).toFixed(1);
    const inv = (1 / transform2.k).toFixed(6);
    return d.symbol === "triangle-down" ? `matrix(-${inv},0,0,-${inv},${tx},${ty})` : `matrix(${inv},0,0,${inv},${tx},${ty})`;
  }
  function connectedMetIdsForRxn(rxnId) {
    const ids = /* @__PURE__ */ new Set();
    if (!fig)
      return ids;
    fig.rail_routes.forEach((rt) => {
      if (rt.rxn_id === rxnId)
        ids.add(rt.met_id);
    });
    fig.inner_edges.forEach((e) => {
      if (e.rxn_id === rxnId)
        ids.add(e.met_id);
    });
    fig.met_nodes.forEach((m) => {
      if (m.consumers.includes(rxnId) || m.producers.includes(rxnId))
        ids.add(m.id);
    });
    return ids;
  }
  function connectedRxnIdsForMet(metId) {
    const ids = /* @__PURE__ */ new Set();
    if (!fig)
      return ids;
    fig.rail_routes.forEach((rt) => {
      if (rt.met_id === metId)
        ids.add(rt.rxn_id);
    });
    fig.inner_edges.forEach((e) => {
      if (e.met_id === metId)
        ids.add(e.rxn_id);
    });
    const met = fig.met_nodes.find((m) => m.id === metId);
    met?.consumers.forEach((r) => ids.add(r));
    met?.producers.forEach((r) => ids.add(r));
    return ids;
  }
  function reactionFocusSet(rxnId) {
    return {
      rxnIds: /* @__PURE__ */ new Set([rxnId]),
      metIds: connectedMetIdsForRxn(rxnId)
    };
  }
  function setRxnFocus(d, ev) {
    ev.stopPropagation();
    if (focusState?.type === "rxn" && focusState.id === d.id) {
      clearFocus();
      return;
    }
    const neighborhood = reactionFocusSet(d.id);
    focusState = {
      type: "rxn",
      id: d.id,
      rxnIds: neighborhood.rxnIds,
      metIds: neighborhood.metIds
    };
    applyFocus();
  }
  function setMetFocus(d, ev) {
    ev.stopPropagation();
    if (focusState?.type === "met" && focusState.id === d.id) {
      clearFocus();
      return;
    }
    focusState = {
      type: "met",
      id: d.id,
      rxnIds: connectedRxnIdsForMet(d.id),
      metIds: /* @__PURE__ */ new Set([d.id])
    };
    applyFocus();
  }
  function focusMetaboliteById(metId) {
    focusState = {
      type: "met",
      id: metId,
      rxnIds: connectedRxnIdsForMet(metId),
      metIds: /* @__PURE__ */ new Set([metId])
    };
    applyFocus();
  }
  function setStationFocus(st, ev) {
    ev.stopPropagation();
    hideTT();
    const rxnIds = /* @__PURE__ */ new Set();
    const metIds = /* @__PURE__ */ new Set();
    st.lines.forEach((ln) => {
      ln.rxn_ids.forEach((rid) => {
        rxnIds.add(rid);
        connectedMetIdsForRxn(rid).forEach((mid) => metIds.add(mid));
      });
    });
    focusState = {
      type: "station",
      id: st.pair_id,
      rxnIds,
      metIds
    };
    applyFocus();
  }
  function setStationReactionFocus(rxnId, ev) {
    ev.stopPropagation();
    focusStationReactionById(rxnId);
  }
  function focusStationReactionById(rxnId) {
    const neighborhood = reactionFocusSet(rxnId);
    focusState = {
      type: "rxn",
      id: rxnId,
      rxnIds: neighborhood.rxnIds,
      metIds: neighborhood.metIds
    };
    applyFocus();
  }
  function clearFocus() {
    if (!focusState)
      return;
    focusState = null;
    hideStationPanel();
    applyFocus();
  }
  function focusHasStation(st) {
    if (!focusState)
      return false;
    return st.lines.some((ln) => ln.rxn_ids.some((rid) => focusState.rxnIds.has(rid)));
  }
  function flowDirection(rxnId, metId, stoich) {
    if (typeof stoich === "number" && Number.isFinite(stoich) && stoich !== 0) {
      return stoich > 0 ? "rxn-to-met" : "met-to-rxn";
    }
    const met = fig?.met_nodes.find((m) => m.id === metId);
    if (met?.producers.includes(rxnId))
      return "rxn-to-met";
    if (met?.consumers.includes(rxnId))
      return "met-to-rxn";
    return null;
  }
  function flowMagnitude(rxnId, metId, stoich) {
    if (typeof stoich === "number" && Number.isFinite(stoich) && stoich !== 0) {
      return Math.max(1, Math.abs(stoich));
    }
    if (!fig)
      return 1;
    let m = 0;
    for (const r of fig.rail_routes) {
      if (r.rxn_id === rxnId && r.met_id === metId)
        m += Math.abs(r.stoich) || 1;
    }
    if (m === 0) {
      for (const e of fig.inner_edges) {
        if (e.rxn_id === rxnId && e.met_id === metId)
          m += 1;
      }
    }
    return Math.max(1, m);
  }
  function applyFlowArrow(sel, isFocus, direction, magnitude = 1) {
    const activeArrow = focusState !== null && isFocus && direction !== null;
    sel.classed("flov-flow-arrow", activeArrow).classed("flov-flow-forward", activeArrow && direction === "rxn-to-met").classed("flov-flow-reverse", activeArrow && direction === "met-to-rxn").style("animation-duration", activeArrow ? `${1 / magnitude}s` : "").attr("marker-end", null).attr("marker-start", null);
  }
  function applyFocus() {
    if (!fig)
      return;
    const zl = wrapper.querySelector(".zoom-layer");
    if (!zl)
      return;
    const active = focusState !== null;
    const dimStroke = lightMode ? "#8c959f" : "#484f58";
    const dimFill = lightMode ? "#d8dee4" : "#30363d";
    const view = fig.views[currentView];
    const railViewById = /* @__PURE__ */ new Map();
    view.rail_routes.forEach((rt) => railViewById.set(rt.line_id, rt));
    select_default2(zl).select(".g-guide").selectAll("path").attr("opacity", (g) => !active || g.node_type === "rxn" && focusState.rxnIds.has(g.node_id) || g.node_type === "met" && focusState.metIds.has(g.node_id) ? 1 : 0.08);
    select_default2(zl).select(".g-inner").selectAll("path").each(function(e) {
      const isFocus = active && focusState.rxnIds.has(e.rxn_id) && focusState.metIds.has(e.met_id);
      const sel = select_default2(this);
      sel.attr("stroke", !active || isFocus ? "rgba(139,148,158,0.50)" : dimStroke).attr("opacity", !active || isFocus ? 1 : 0.1);
      applyFlowArrow(sel, isFocus, flowDirection(e.rxn_id, e.met_id), flowMagnitude(e.rxn_id, e.met_id));
    });
    select_default2(zl).select(".g-rail").selectAll("path.rail-line").each(function(d) {
      const rv = railViewById.get(d.line_id);
      const isFocus = active && focusState.rxnIds.has(d.rxn_id) && focusState.metIds.has(d.met_id);
      const sel = select_default2(this);
      sel.attr("stroke", !active || isFocus ? rv?.color ?? "rgba(110,118,129,0.18)" : dimStroke).attr("stroke-width", rv?.width ?? 0.8).attr("opacity", !active || isFocus ? 1 : 0.12);
      applyFlowArrow(sel, isFocus, flowDirection(d.rxn_id, d.met_id, d.stoich), flowMagnitude(d.rxn_id, d.met_id, d.stoich));
    });
    select_default2(zl).select(".g-route").selectAll("path.route-line").each(function(d) {
      const sv = view.stations.find((s) => s.pair_id === d.pair_id);
      const lv = sv ? sv.lines[d.lineIndex] : void 0;
      const isFocus = active && (d.rxn_ids ?? []).some((rid) => focusState.rxnIds.has(rid));
      select_default2(this).classed("flov-flow-arrow", false).classed("flov-flow-forward", false).classed("flov-flow-reverse", false).style("animation-duration", "").attr("marker-end", null).attr("marker-start", null).attr("stroke", !active || isFocus ? lv?.color ?? "#6e7681" : dimStroke).attr("stroke-width", lv?.width ?? 1).attr("opacity", !active || isFocus ? 1 : 0.1);
    });
    select_default2(zl).select(".g-stn").selectAll("rect.stn-pill").each(function(p) {
      const stroke = p.side === "a" ? p.st.comp_a_color : p.st.comp_b_color;
      const isFocus = active && focusHasStation(p.st);
      select_default2(this).attr("fill", !active || isFocus ? "#ffffff" : dimFill).attr("stroke", !active || isFocus ? stroke : dimStroke).attr("opacity", !active || isFocus ? 1 : 0.16);
    });
    gMetSel?.selectAll("path.met-shape").each(function(d) {
      const isFocus = active && focusState.metIds.has(d.id);
      select_default2(this).attr("fill", !active || isFocus ? d.display_kind === "stop" ? "#ffffff" : d.comp_color : dimFill).attr("fill-opacity", !active || isFocus ? d.display_kind === "stop" ? 1 : 0.85 : 0.55).attr("stroke", !active || isFocus ? d.display_kind === "stop" ? d.comp_color : "#0d1117" : dimStroke).attr("opacity", !active || isFocus ? 1 : 0.18);
    });
    gRxnSel?.each(function(d) {
      const isFocus = active && focusState.rxnIds.has(d.id);
      select_default2(this).attr("fill", !active || isFocus ? d.fill_color : dimFill).attr("stroke", !active || isFocus ? d.border_color : dimStroke).attr("stroke-width", d.border_width).attr("opacity", !active || isFocus ? d.opacity : 0.16);
    });
    select_default2(zl).select(".g-lbl").selectAll("text").attr("opacity", (d) => !active || focusState.rxnIds.has(d.id) ? 1 : 0.1);
  }
  function buildToolbar(data, parent) {
    const bar = document.createElement("div");
    bar.className = "flov-toolbar";
    parent.appendChild(bar);
    if (data.view_labels.length) {
      const grp = btnGroup(bar, "Flux view");
      data.view_labels.forEach((lbl, i) => {
        const btn = mkBtn(lbl, i === 0);
        btn.onclick = () => switchView(i);
        grp.appendChild(btn);
        viewBtns.push(btn);
      });
    }
    mkToggle(bar, "Metabolites", "Show", "Hide", true, (on) => {
      gMetSel.style("display", on ? "" : "none");
    });
    const themeGroup = btnGroup(bar, "Theme");
    const themeBtn = mkBtn("Dark mode", true);
    const applyTheme = () => {
      wrapper.classList.toggle("light-mode", lightMode);
      if (svgEl)
        svgEl.style.background = lightMode ? "#f6f8fa" : "#0d1117";
      themeBtn.textContent = lightMode ? "Dark mode" : "Light mode";
      themeBtn.classList.toggle("active", lightMode);
      applyFocus();
    };
    themeBtn.onclick = () => {
      lightMode = !lightMode;
      applyTheme();
    };
    applyTheme();
    themeGroup.appendChild(themeBtn);
    const rg = btnGroup(bar, "");
    const rb = mkBtn("\u2316 Reset zoom", false);
    rb.onclick = () => {
      if (svgEl && zoomB)
        select_default2(svgEl).transition().duration(350).call(zoomB.transform, identity3);
    };
    rg.appendChild(rb);
  }
  function btnGroup(bar, label) {
    const grp = document.createElement("div");
    grp.className = "flov-btn-group";
    if (label) {
      const lbl = document.createElement("span");
      lbl.className = "flov-btn-label";
      lbl.textContent = label;
      grp.appendChild(lbl);
    }
    bar.appendChild(grp);
    return grp;
  }
  function mkBtn(text, active) {
    const b = document.createElement("button");
    b.className = "flov-btn" + (active ? " active" : "");
    b.textContent = text;
    return b;
  }
  function mkToggle(bar, label, onTxt, offTxt, initOn, cb) {
    const grp = btnGroup(bar, label);
    const bOn = mkBtn(onTxt, initOn);
    const bOff = mkBtn(offTxt, !initOn);
    bOn.onclick = () => {
      bOn.classList.add("active");
      bOff.classList.remove("active");
      cb(true);
    };
    bOff.onclick = () => {
      bOff.classList.add("active");
      bOn.classList.remove("active");
      cb(false);
    };
    grp.appendChild(bOn);
    grp.appendChild(bOff);
  }
  function switchView(vi) {
    if (!fig)
      return;
    currentView = vi;
    const view = fig.views[vi];
    const zl = wrapper.querySelector(".zoom-layer");
    if (!zl)
      return;
    const railViewById = /* @__PURE__ */ new Map();
    view.rail_routes.forEach((rt) => railViewById.set(rt.line_id, rt));
    select_default2(zl).select(".g-rail").selectAll("path.rail-line").each(function(d) {
      const rv = railViewById.get(d.line_id);
      select_default2(this).attr("stroke", rv?.color ?? "rgba(110,118,129,0.18)").attr("stroke-width", rv?.width ?? 0.8);
    });
    select_default2(zl).select(".g-rxn").selectAll("path.rxn-shape").data(view.rxn_nodes).attr("d", (d) => symPath(d.symbol, d.r)).attr("transform", (d) => rxnTransform(d)).attr("fill", (d) => d.fill_color).attr("opacity", (d) => d.opacity);
    select_default2(zl).select(".g-route").selectAll("path.route-line").each(function(d) {
      const sv = view.stations.find((s) => s.pair_id === d.pair_id);
      const lv = sv ? sv.lines[d.lineIndex] : void 0;
      if (lv) {
        select_default2(this).attr("d", railOffsetPath(d.path, d.seg_slots, xS, yS, currentTransform.k)).attr("stroke", lv.color).attr("stroke-width", lv.width);
      } else {
        select_default2(this).attr("d", railOffsetPath(d.path, d.seg_slots, xS, yS, currentTransform.k)).attr("stroke", "#6e7681").attr("stroke-width", 1);
      }
    });
    gRxnHitSel.data(view.rxn_nodes).attr("d", (d) => symPath(d.symbol, d.r)).attr("transform", (d) => rxnTransform(d));
    select_default2(zl).select(".g-lbl").selectAll("text").data(view.rxn_nodes).attr("x", (d) => xS(d.x)).attr("y", (d) => yS(d.y) - d.r - 2);
    viewBtns.forEach((b, i) => b.classList.toggle("active", i === vi));
    applyFocus();
    refreshStationPanelForView();
  }
  function stationClassLabel(klass) {
    return klass.replace(/_/g, " ").replace(/\b\w/g, (c) => c.toUpperCase());
  }
  function compactMetList(mets) {
    if (!mets || !mets.length)
      return "---";
    const text = mets.slice(0, 2).join(", ");
    return mets.length > 2 ? `${text}, ...` : text;
  }
  function activeFluxLabel(rxn) {
    const raw = rxn.flux_per_view[currentView] ?? rxn.flux_per_view[0] ?? "";
    const idx = raw.indexOf(":");
    return idx >= 0 ? raw.slice(idx + 1).trim() : raw;
  }
  function stationReactionSearchText(rxn) {
    return [
      rxn.id,
      rxn.name,
      ...rxn.substrates ?? [],
      ...rxn.products ?? [],
      ...rxn.flux_per_view ?? []
    ].join(" ").toLowerCase();
  }
  function buildStationPanel(parent) {
    const el = document.createElement("div");
    el.className = "flov-station-panel";
    el.addEventListener("click", (ev) => ev.stopPropagation());
    parent.appendChild(el);
    stationPanelEl = el;
  }
  function hideStationPanel() {
    if (!stationPanelEl)
      return;
    stationPanelEl.style.display = "none";
    stationPanelEl.innerHTML = "";
    activeStationPanel = null;
  }
  function refreshStationPanelForView() {
    if (!stationPanelEl || stationPanelEl.style.display === "none" || !activeStationPanel) {
      return;
    }
    showStationPanel(activeStationPanel, true);
  }
  function showStationPanel(st, preserveState = false) {
    if (!stationPanelEl)
      return;
    const scrollTop = preserveState ? stationPanelEl.scrollTop : 0;
    const openGroups = preserveState ? new Set(
      Array.from(stationPanelEl.querySelectorAll(".flov-station-group.open")).map((group) => group.dataset.stationClass ?? "").filter(Boolean)
    ) : /* @__PURE__ */ new Set();
    activeStationPanel = st;
    stationPanelEl.innerHTML = "";
    stationPanelEl.style.borderColor = st.comp_a_color;
    const head = document.createElement("div");
    head.className = "flov-station-head";
    const title = document.createElement("div");
    title.className = "flov-station-title";
    const titleText = document.createElement("span");
    titleText.textContent = `${st.comp_a_label} \u2194 ${st.comp_b_label}`;
    const count = document.createElement("span");
    count.className = "flov-station-badge";
    count.textContent = `${st.n_rxns} rxn${st.n_rxns === 1 ? "" : "s"}`;
    const close = document.createElement("button");
    close.className = "flov-station-close";
    close.type = "button";
    close.setAttribute("aria-label", "Close station flows");
    close.textContent = "\xD7";
    close.onclick = (ev) => {
      ev.stopPropagation();
      if (focusState?.type === "station" && focusState.id === st.pair_id) {
        clearFocus();
      } else {
        hideStationPanel();
      }
    };
    title.appendChild(titleText);
    title.appendChild(count);
    title.appendChild(close);
    const sub = document.createElement("div");
    sub.className = "flov-station-sub";
    sub.textContent = `Reaction flows \xB7 ${fig?.views[currentView]?.label ?? ""}`;
    head.appendChild(title);
    head.appendChild(sub);
    const body = document.createElement("div");
    body.className = "flov-station-body";
    st.lines.forEach((ln) => {
      const group = document.createElement("div");
      group.className = "flov-station-group";
      group.dataset.stationClass = ln.klass;
      const toggle = document.createElement("button");
      toggle.className = "flov-station-toggle";
      toggle.type = "button";
      const caret = document.createElement("span");
      caret.className = "flov-station-caret";
      caret.textContent = "\u25B8";
      const swatch = document.createElement("span");
      swatch.className = "flov-station-swatch";
      swatch.style.background = KLASS_COLORS[ln.klass] ?? "#6e7681";
      const label = document.createElement("span");
      label.textContent = stationClassLabel(ln.klass);
      const rxnCount = document.createElement("span");
      rxnCount.className = "flov-station-count";
      rxnCount.textContent = `${ln.rxn_summaries.length}`;
      toggle.appendChild(caret);
      toggle.appendChild(swatch);
      toggle.appendChild(label);
      toggle.appendChild(rxnCount);
      const rxns = document.createElement("div");
      rxns.className = "flov-station-rxns";
      appendStationRailMaps(rxns, ln, st);
      toggle.onclick = (ev) => {
        ev.stopPropagation();
        const open = group.classList.toggle("open");
        caret.textContent = open ? "\u25BE" : "\u25B8";
      };
      group.appendChild(toggle);
      group.appendChild(rxns);
      body.appendChild(group);
    });
    stationPanelEl.appendChild(head);
    stationPanelEl.appendChild(body);
    stationPanelEl.style.display = "block";
    if (preserveState) {
      stationPanelEl.querySelectorAll(".flov-station-group").forEach((group) => {
        if (!openGroups.has(group.dataset.stationClass ?? ""))
          return;
        group.classList.add("open");
        const caret = group.querySelector(".flov-station-caret");
        if (caret)
          caret.textContent = "\u25BE";
      });
      stationPanelEl.scrollTop = scrollTop;
    } else {
      stationPanelEl.scrollTop = 0;
    }
    applyStationPanelSearch(currentStationSearch);
  }
  function applyStationPanelSearch(query, scrollToRxnId) {
    if (!stationPanelEl || stationPanelEl.style.display === "none")
      return;
    const lq = query.trim().toLowerCase();
    const shouldFilter = lq.length > 0;
    stationPanelEl.querySelectorAll(".flov-station-group").forEach((group) => {
      const branches = Array.from(group.querySelectorAll(".flov-station-branch"));
      let groupHasMatch = false;
      branches.forEach((branch) => {
        const isTarget = scrollToRxnId !== void 0 && branch.dataset.rxnId === scrollToRxnId;
        const isMatch = !shouldFilter || (branch.dataset.searchText ?? "").includes(lq);
        groupHasMatch ||= isMatch || isTarget;
        branch.classList.toggle("search-hidden", shouldFilter && !isMatch && !isTarget);
        branch.classList.toggle("search-hit", shouldFilter && isMatch || isTarget);
      });
      group.querySelectorAll(".flov-station-map").forEach((svg) => {
        const hasVisibleBranch = Array.from(
          svg.querySelectorAll(".flov-station-branch")
        ).some((branch) => !branch.classList.contains("search-hidden"));
        svg.classList.toggle("search-hidden", shouldFilter && !hasVisibleBranch);
      });
      group.classList.toggle("search-hidden", shouldFilter && !groupHasMatch);
      const caret = group.querySelector(".flov-station-caret");
      if (shouldFilter && groupHasMatch) {
        group.classList.add("open");
        if (caret)
          caret.textContent = "\u25BE";
      } else if (!shouldFilter) {
        branches.forEach((branch) => branch.classList.remove("search-hit", "search-hidden"));
        group.querySelectorAll(".flov-station-map").forEach((svg) => {
          svg.classList.remove("search-hidden");
        });
        group.classList.remove("search-hidden");
      }
    });
    if (scrollToRxnId) {
      const target = Array.from(
        stationPanelEl.querySelectorAll(".flov-station-branch")
      ).find((branch) => branch.dataset.rxnId === scrollToRxnId);
      target?.scrollIntoView({ block: "center", inline: "nearest" });
    }
  }
  function appendStationRailMaps(parent, line, station) {
    const svgNS = "http://www.w3.org/2000/svg";
    const color2 = KLASS_COLORS[line.klass] ?? "#6e7681";
    const chunkSize = 5;
    for (let start2 = 0; start2 < line.rxn_summaries.length; start2 += chunkSize) {
      const chunk = line.rxn_summaries.slice(start2, start2 + chunkSize);
      const width = 650;
      const height = 225;
      const railY = 160;
      const x0 = 112;
      const splitStartX = 205;
      const splitGap = 62;
      const x1 = width - 126;
      const laneGap = 9;
      const markerPad = 4;
      const laneTop = railY - (chunk.length - 1) / 2 * laneGap;
      const laneBottom = railY + (chunk.length - 1) / 2 * laneGap;
      const bundleTop = laneTop - markerPad;
      const bundleBottom = laneBottom + markerPad;
      const svg = document.createElementNS(svgNS, "svg");
      svg.setAttribute("class", "flov-station-map");
      svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
      svg.setAttribute("role", "img");
      svg.setAttribute(
        "aria-label",
        `${stationClassLabel(line.klass)} reaction rail ${start2 / chunkSize + 1}`
      );
      const endLaneY = railY + (chunk.length - 1 - (chunk.length - 1) / 2) * laneGap;
      [
        [x0, station.comp_a_color, bundleTop, bundleBottom],
        [x1, station.comp_b_color, endLaneY - markerPad, endLaneY + markerPad]
      ].forEach(([x2, stroke, yTop, yBottom]) => {
        const pill = document.createElementNS(svgNS, "rect");
        pill.setAttribute("x", String(Number(x2) - 4));
        pill.setAttribute("y", String(yTop));
        pill.setAttribute("width", "8");
        pill.setAttribute("height", String(Number(yBottom) - Number(yTop)));
        pill.setAttribute("rx", "4");
        pill.setAttribute("fill", "#fff");
        pill.setAttribute("stroke", String(stroke));
        pill.setAttribute("stroke-width", "2");
        svg.appendChild(pill);
      });
      chunk.forEach((rxn, i) => {
        const isThroughLine = i === chunk.length - 1;
        const laneY = railY + (i - (chunk.length - 1) / 2) * laneGap;
        const splitX = splitStartX + i * splitGap;
        const branchEndX = isThroughLine ? x1 - 34 : Math.min(splitX + 76, x1 - 42);
        const branchEndY = isThroughLine ? laneY : 34 + i * 31;
        const labelX = isThroughLine ? (splitStartX + Math.max(0, chunk.length - 2) * splitGap + x1) / 2 + 42 : (branchEndX + x1) / 2;
        const labelY = branchEndY - (isThroughLine ? 28 : 14);
        const fluxX = isThroughLine ? (splitStartX + Math.max(0, chunk.length - 2) * splitGap + x1) / 2 : ((i === 0 ? x0 : splitStartX + (i - 1) * splitGap) + splitX) / 2;
        const fluxY = laneY - 5;
        const labelAnchor = branchEndX > width * 0.76 ? "end" : "middle";
        const branchColor = STATION_BRANCH_COLORS[i % STATION_BRANCH_COLORS.length];
        const substrate = compactMetList(rxn.substrates);
        const product = compactMetList(rxn.products);
        const g = document.createElementNS(svgNS, "g");
        g.setAttribute("class", "flov-station-branch");
        g.dataset.rxnId = rxn.id;
        g.dataset.searchText = stationReactionSearchText(rxn);
        g.addEventListener("click", (ev) => setStationReactionFocus(rxn.id, ev));
        const subLabel = document.createElementNS(svgNS, "text");
        subLabel.setAttribute("class", "flov-station-map-met");
        subLabel.setAttribute("x", String(x0 - 10));
        subLabel.setAttribute("y", (laneY + 3).toFixed(1));
        subLabel.setAttribute("text-anchor", "end");
        subLabel.textContent = substrate.length > 20 ? `${substrate.slice(0, 19)}...` : substrate;
        g.appendChild(subLabel);
        const branch = document.createElementNS(svgNS, "path");
        branch.setAttribute(
          "d",
          isThroughLine ? `M${x0} ${laneY.toFixed(1)} L${x1} ${laneY.toFixed(1)}` : `M${x0} ${laneY.toFixed(1)} L${splitX} ${laneY.toFixed(1)} L${branchEndX.toFixed(1)} ${branchEndY.toFixed(1)} L${x1} ${branchEndY.toFixed(1)}`
        );
        branch.setAttribute("class", "flov-station-map-line");
        branch.setAttribute("stroke", branchColor);
        branch.setAttribute("stroke-width", "1.9");
        branch.setAttribute("fill", "none");
        g.appendChild(branch);
        const productLabel = document.createElementNS(svgNS, "text");
        productLabel.setAttribute("class", "flov-station-map-met");
        productLabel.setAttribute("x", String(x1 + 10));
        productLabel.setAttribute("y", (branchEndY + 3).toFixed(1));
        productLabel.setAttribute("text-anchor", "start");
        productLabel.textContent = product.length > 20 ? `${product.slice(0, 19)}...` : product;
        g.appendChild(productLabel);
        if (!isThroughLine) {
          const dot = document.createElementNS(svgNS, "circle");
          dot.setAttribute("cx", String(splitX));
          dot.setAttribute("cy", laneY.toFixed(1));
          dot.setAttribute("r", "2.7");
          dot.setAttribute("fill", "#fff");
          dot.setAttribute("stroke", branchColor);
          dot.setAttribute("stroke-width", "1.5");
          g.appendChild(dot);
        }
        const label = document.createElementNS(svgNS, "text");
        label.setAttribute("class", "flov-station-map-label");
        label.setAttribute("x", labelX.toFixed(1));
        label.setAttribute("y", String(labelY));
        label.setAttribute("text-anchor", labelAnchor);
        label.textContent = rxn.id.length > 18 ? `${rxn.id.slice(0, 17)}...` : rxn.id;
        g.appendChild(label);
        const flux = document.createElementNS(svgNS, "text");
        flux.setAttribute("class", "flov-station-map-flux");
        flux.setAttribute("x", fluxX.toFixed(1));
        flux.setAttribute("y", String(fluxY));
        flux.setAttribute("text-anchor", "middle");
        flux.textContent = activeFluxLabel(rxn);
        g.appendChild(flux);
        const title = document.createElementNS(svgNS, "title");
        title.textContent = `${rxn.id}${rxn.name && rxn.name !== rxn.id ? ` \xB7 ${rxn.name}` : ""}
Sub: ${compactMetList(rxn.substrates)}
Prd: ${compactMetList(rxn.products)}
${rxn.flux_per_view.join("\n")}`;
        g.appendChild(title);
        svg.appendChild(g);
      });
      chunk.slice(0, -1).forEach((_, i) => {
        const splitX = splitStartX + i * splitGap;
        const remainingLaneYs = chunk.slice(i).map((__, ri) => railY + (i + ri - (chunk.length - 1) / 2) * laneGap);
        const markerTop = Math.min(...remainingLaneYs) - markerPad;
        const markerBottom = Math.max(...remainingLaneYs) + markerPad;
        const marker = document.createElementNS(svgNS, "rect");
        marker.setAttribute("x", String(splitX - 3.5));
        marker.setAttribute("y", String(markerTop));
        marker.setAttribute("width", "7");
        marker.setAttribute("height", String(markerBottom - markerTop));
        marker.setAttribute("rx", "3.5");
        marker.setAttribute("fill", "#fff");
        marker.setAttribute("stroke", color2);
        marker.setAttribute("stroke-width", "1.5");
        svg.appendChild(marker);
      });
      parent.appendChild(svg);
    }
  }
  function buildTooltip(parent) {
    const el = document.createElement("div");
    el.className = "flov-tt";
    parent.appendChild(el);
    tooltipEl = el;
  }
  function showTT(ev, html, borderColor) {
    if (!tooltipEl)
      return;
    tooltipEl.innerHTML = html;
    tooltipEl.style.borderColor = borderColor;
    tooltipEl.style.display = "block";
    moveTT(ev);
  }
  function moveTT(ev) {
    if (!tooltipEl)
      return;
    const rect = wrapper.getBoundingClientRect();
    let tx = ev.clientX - rect.left + 14;
    let ty = ev.clientY - rect.top - 10;
    const tw = tooltipEl.offsetWidth || 220;
    const th = tooltipEl.offsetHeight || 110;
    if (tx + tw > rect.width - 6)
      tx = ev.clientX - rect.left - tw - 14;
    if (ty + th > rect.height - 6)
      ty = ev.clientY - rect.top - th - 10;
    tooltipEl.style.left = tx + "px";
    tooltipEl.style.top = ty + "px";
  }
  function hideTT() {
    if (tooltipEl)
      tooltipEl.style.display = "none";
  }
  function buildLegend(parent, data, svgH) {
    const SVG_NS = "http://www.w3.org/2000/svg";
    const gradId = `flov-legend-grad-${Math.random().toString(36).slice(2)}`;
    const panel = document.createElement("div");
    panel.className = "flov-legend";
    panel.style.maxHeight = Math.max(svgH - 70, 280) + "px";
    parent.appendChild(panel);
    const head = document.createElement("div");
    head.className = "flov-legend-head";
    const headTitle = document.createElement("span");
    headTitle.className = "flov-legend-title";
    headTitle.textContent = "Legend";
    const headCaret = document.createElement("span");
    headCaret.className = "flov-legend-caret";
    headCaret.textContent = "\u25BE";
    head.appendChild(headTitle);
    head.appendChild(headCaret);
    panel.appendChild(head);
    const body = document.createElement("div");
    body.className = "flov-legend-body";
    panel.appendChild(body);
    head.onclick = () => {
      const open = !panel.classList.contains("collapsed");
      panel.classList.toggle("collapsed", open);
      headCaret.textContent = open ? "\u25B8" : "\u25BE";
    };
    function section(title) {
      const sec = document.createElement("div");
      sec.className = "flov-legend-section";
      const t = document.createElement("div");
      t.className = "flov-legend-section-title";
      t.textContent = title;
      sec.appendChild(t);
      body.appendChild(sec);
      return sec;
    }
    function row(sec, swatch, label, sub) {
      const r = document.createElement("div");
      r.className = "flov-legend-row";
      const sw = document.createElement("div");
      sw.className = "flov-legend-swatch";
      sw.appendChild(swatch);
      const text = document.createElement("div");
      text.className = "flov-legend-text";
      const main = document.createElement("div");
      main.className = "flov-legend-label";
      main.textContent = label;
      text.appendChild(main);
      if (sub) {
        const s = document.createElement("div");
        s.className = "flov-legend-sub";
        s.textContent = sub;
        text.appendChild(s);
      }
      r.appendChild(sw);
      r.appendChild(text);
      sec.appendChild(r);
    }
    function mkSvg(w, h) {
      const s = document.createElementNS(SVG_NS, "svg");
      s.setAttribute("width", String(w));
      s.setAttribute("height", String(h));
      s.setAttribute("viewBox", `0 0 ${w} ${h}`);
      return s;
    }
    function mkEl(tag, attrs) {
      const e = document.createElementNS(SVG_NS, tag);
      for (const [k, v] of Object.entries(attrs))
        e.setAttribute(k, String(v));
      return e;
    }
    {
      const sec = section("Flux (mmol/gDW/h)");
      const cbW = 198, cbH = 38;
      const svg = mkSvg(cbW, cbH);
      const defs = mkEl("defs", {});
      const grad = mkEl("linearGradient", { id: gradId, x1: "0", y1: "0", x2: "1", y2: "0" });
      const stops = [
        ["0%", "#1565c0"],
        ["25%", "#90caf9"],
        ["50%", "#6e7681"],
        ["75%", "#ef9a9a"],
        ["100%", "#c62828"]
      ];
      stops.forEach(([off, col]) => {
        grad.appendChild(mkEl("stop", { offset: off, "stop-color": col }));
      });
      defs.appendChild(grad);
      svg.appendChild(defs);
      svg.appendChild(mkEl("rect", { x: 2, y: 4, width: cbW - 4, height: 12, rx: 3, fill: `url(#${gradId})` }));
      const absMax = data.meta.abs_max_flux;
      const fmt = format("+.2f");
      const labels = [
        [2, fmt(-absMax)],
        [cbW / 2, "0"],
        [cbW - 2, fmt(absMax)]
      ];
      labels.forEach(([x2, txt], i) => {
        const t = mkEl("text", {
          x: x2,
          y: 28,
          "text-anchor": i === 0 ? "start" : i === labels.length - 1 ? "end" : "middle",
          "font-size": 9,
          "font-family": "'Inter',Arial,sans-serif",
          fill: "#8b949e"
        });
        t.textContent = txt;
        svg.appendChild(t);
      });
      const wrap = document.createElement("div");
      wrap.className = "flov-legend-fluxbar";
      wrap.appendChild(svg);
      sec.appendChild(wrap);
      const note = document.createElement("div");
      note.className = "flov-legend-note";
      note.textContent = "Grey = inactive / no data";
      sec.appendChild(note);
    }
    {
      const sec = section("Reactions");
      const sv1 = mkSvg(28, 20);
      sv1.appendChild(mkEl("circle", { cx: 14, cy: 10, r: 7, fill: "#ef9a9a", stroke: "#58a6ff", "stroke-width": 2 }));
      row(sec, sv1, "Reaction node", "Fill = flux  \xB7  Border = compartment");
      const sv2 = mkSvg(28, 20);
      sv2.appendChild(mkEl("circle", { cx: 14, cy: 10, r: 7, fill: "#6e7681", stroke: "#58a6ff", "stroke-width": 2, opacity: 0.22 }));
      row(sec, sv2, "No flux data with this solver", "Faded grey fill");
    }
    {
      const sec = section("Metabolites");
      const sv1 = mkSvg(28, 20);
      sv1.appendChild(mkEl("circle", { cx: 14, cy: 10, r: 6.5, fill: "#ffffff", stroke: "#58a6ff", "stroke-width": 1.6 }));
      row(sec, sv1, "Metabolite", "Stop on a rail line");
      const sv2 = mkSvg(28, 20);
      sv2.appendChild(mkEl("path", {
        d: "M14,3 L21,10 L14,17 L7,10 Z",
        fill: "#58a6ff",
        "fill-opacity": 0.85,
        stroke: "#0d1117",
        "stroke-width": 0.7
      }));
      row(sec, sv2, "Detached metabolite", "Only connected to other compartment");
    }
    if (data.compartments && data.compartments.length) {
      const sec = section("Compartments");
      for (const c of data.compartments) {
        const sv = mkSvg(28, 20);
        sv.appendChild(mkEl("rect", {
          x: 3,
          y: 4,
          width: 22,
          height: 12,
          rx: 3,
          fill: c.fill,
          stroke: c.color,
          "stroke-width": 1.4,
          "stroke-dasharray": "4,2"
        }));
        row(sec, sv, c.label || c.key);
      }
    }
    {
      const sec = section("Stations");
      const sv = mkSvg(32, 18);
      sv.appendChild(mkEl("rect", {
        x: 4,
        y: 5,
        width: 24,
        height: 8,
        rx: 4,
        ry: 4,
        fill: "#ffffff",
        stroke: "#58a6ff",
        "stroke-width": 2.2
      }));
      row(sec, sv, "Station hub", "Compartment interconnect");
    }
    {
      const sec = section("Station line classes");
      const entries = [
        ["amino_acid", "Amino acids"],
        ["sugar", "Sugars"],
        ["cofactor", "Cofactors"],
        ["inorganic", "Inorganics"],
        ["other", "Other"]
      ];
      for (const [k, lbl] of entries) {
        const sv = mkSvg(28, 14);
        sv.appendChild(mkEl("line", {
          x1: 2,
          y1: 7,
          x2: 26,
          y2: 7,
          stroke: KLASS_COLORS[k] ?? "#ec4899",
          "stroke-width": 3.4,
          "stroke-linecap": "round"
        }));
        row(sec, sv, lbl);
      }
      const svI = mkSvg(28, 14);
      svI.appendChild(mkEl("line", {
        x1: 2,
        y1: 7,
        x2: 26,
        y2: 7,
        stroke: "#6e7681",
        "stroke-width": 2.2,
        "stroke-linecap": "round"
      }));
      row(sec, svI, "Inactive class", "No flux on this line");
    }
    {
      const sec = section("Lines");
      const sv1 = mkSvg(28, 14);
      sv1.appendChild(mkEl("line", {
        x1: 2,
        y1: 7,
        x2: 26,
        y2: 7,
        stroke: "#58a6ff",
        "stroke-width": 2.6,
        "stroke-linecap": "round"
      }));
      row(sec, sv1, "Active line", "Solid: flux on this edge");
      const sv2 = mkSvg(28, 14);
      const ln2 = mkEl("line", {
        x1: 2,
        y1: 7,
        x2: 26,
        y2: 7,
        stroke: "#58a6ff",
        "stroke-width": 2.6,
        "stroke-linecap": "round"
      });
      ln2.setAttribute("class", "flov-flow-arrow flov-flow-forward");
      sv2.appendChild(ln2);
      row(sec, sv2, "Flow: rxn \u2192 met", "Produced by reaction");
      const sv3 = mkSvg(28, 14);
      const ln3 = mkEl("line", {
        x1: 2,
        y1: 7,
        x2: 26,
        y2: 7,
        stroke: "#58a6ff",
        "stroke-width": 2.6,
        "stroke-linecap": "round"
      });
      ln3.setAttribute("class", "flov-flow-arrow flov-flow-reverse");
      sv3.appendChild(ln3);
      row(sec, sv3, "Flow: met \u2192 rxn", "Consumed by reaction");
      const sv4 = mkSvg(28, 14);
      sv4.appendChild(mkEl("line", {
        x1: 2,
        y1: 7,
        x2: 26,
        y2: 7,
        stroke: "#6e7681",
        "stroke-width": 2.2,
        "stroke-linecap": "round",
        opacity: 0.6
      }));
      row(sec, sv4, "Inactive", "No flux with this solver");
    }
  }
  function buildSearch(parent, data) {
    const ia = data.interactivity;
    const rxnIdsLc = ia.rxn_ids.map((s) => s.toLowerCase());
    const rxnNamesLc = ia.rxn_names.map((s) => s.toLowerCase());
    const metIdsLc = ia.met_ids.map((s) => s.toLowerCase());
    const metNamesLc = ia.met_names.map((s) => s.toLowerCase());
    const wrap = document.createElement("div");
    wrap.className = "flov-search-wrap";
    const inp = document.createElement("input");
    inp.className = "flov-search-input";
    inp.placeholder = "\u{1F50D} Search reaction / metabolite\u2026";
    inp.type = "text";
    const res = document.createElement("div");
    res.className = "flov-search-results";
    wrap.appendChild(inp);
    wrap.appendChild(res);
    parent.appendChild(wrap);
    function selectResult(type2, idx) {
      inp.value = type2 === "rxn" ? ia.rxn_names[idx] || ia.rxn_ids[idx] : ia.met_names[idx] || ia.met_ids[idx];
      res.style.display = "none";
      if (type2 === "rxn") {
        focusStationReactionById(ia.rxn_ids[idx]);
      } else {
        focusMetaboliteById(ia.met_ids[idx]);
      }
    }
    function appendSearchResultText(el, id2, name, suffix) {
      const idEl = document.createElement("b");
      idEl.textContent = id2;
      el.appendChild(idEl);
      if (name && name !== id2) {
        const nameEl = document.createElement("span");
        nameEl.style.color = "#7d8590";
        nameEl.textContent = ` \u2014 ${name}`;
        el.appendChild(nameEl);
      }
      const suffixEl = document.createElement("span");
      suffixEl.style.color = "#444c56";
      suffixEl.style.fontSize = "10px";
      suffixEl.textContent = ` ${suffix}`;
      el.appendChild(suffixEl);
    }
    function selectStationResult(hit) {
      showStationPanel(hit.station);
      focusStationReactionById(hit.rxn.id);
      inp.value = hit.rxn.name && hit.rxn.name !== hit.rxn.id ? hit.rxn.name : hit.rxn.id;
      currentStationSearch = inp.value;
      applyStationPanelSearch(currentStationSearch, hit.rxn.id);
      res.style.display = "none";
    }
    function runSearch(q) {
      res.innerHTML = "";
      currentStationSearch = q;
      applyStationPanelSearch(q);
      if (!q.trim()) {
        res.style.display = "none";
        return;
      }
      const lq = q.toLowerCase();
      const hits = [];
      for (let i = 0; i < ia.rxn_ids.length; i++) {
        if (rxnIdsLc[i].includes(lq) || rxnNamesLc[i].includes(lq)) {
          if (hits.length < 20)
            hits.push({ type: "rxn", idx: i });
        }
      }
      for (let i = 0; i < ia.met_ids.length; i++) {
        if (metIdsLc[i].includes(lq) || metNamesLc[i].includes(lq)) {
          if (hits.length < 20)
            hits.push({ type: "met", idx: i });
        }
      }
      for (const station of data.stations) {
        for (const line of station.lines) {
          for (const rxn of line.rxn_summaries) {
            if (!stationReactionSearchText(rxn).includes(lq))
              continue;
            if (hits.length < 20)
              hits.push({ type: "station-rxn", station, line, rxn });
          }
        }
      }
      if (!hits.length) {
        const d = document.createElement("div");
        d.className = "flov-search-result";
        d.style.color = "#7d8590";
        d.textContent = "No matches";
        res.appendChild(d);
        res.style.display = "block";
        return;
      }
      for (const h of hits) {
        const el = document.createElement("div");
        el.className = "flov-search-result";
        if (h.type === "station-rxn") {
          const id2 = h.rxn.id;
          const name = h.rxn.name;
          el.title = `${id2}${name && name !== id2 ? " \u2014 " + name : ""}`;
          appendSearchResultText(el, id2, name, `[station ${stationClassLabel(h.line.klass)}]`);
          el.onclick = () => selectStationResult(h);
        } else {
          const id2 = h.type === "rxn" ? ia.rxn_ids[h.idx] : ia.met_ids[h.idx];
          const name = h.type === "rxn" ? ia.rxn_names[h.idx] : ia.met_names[h.idx];
          el.title = id2 + (name && name !== id2 ? " \u2014 " + name : "");
          appendSearchResultText(el, id2, name, `[${h.type}]`);
          el.onclick = () => selectResult(h.type, h.idx);
        }
        res.appendChild(el);
      }
      res.style.display = "block";
    }
    inp.addEventListener("input", function() {
      runSearch(this.value);
    });
    inp.addEventListener("keydown", (e) => {
      if (e.key === "Escape") {
        inp.value = "";
        currentStationSearch = "";
        applyStationPanelSearch("");
        if (focusState) {
          focusState = null;
          applyFocus();
        }
        res.style.display = "none";
      }
    });
    document.addEventListener("click", function outside(e) {
      if (!wrap.isConnected) {
        document.removeEventListener("click", outside);
        return;
      }
      if (!wrap.contains(e.target))
        res.style.display = "none";
    });
  }
  new ResizeObserver((entries) => {
    for (const entry of entries) {
      const w = entry.contentRect.width;
      if (w <= 0)
        continue;
      if (!fig)
        continue;
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(() => buildUI(fig), 120);
    }
  }).observe(I);
  function onData() {
    const raw = S.get("_figure_json");
    if (!raw || raw === "{}")
      return;
    let data;
    try {
      data = JSON.parse(raw);
    } catch {
      return;
    }
    buildUI(data);
  }
  S.on("change:_figure_json", onData);
  onData();
}
var flux_network_default = { render };
export {
  flux_network_default as default
};
