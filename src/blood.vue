<script setup>
import {nextTick, onMounted, reactive, ref, toRef, toRefs} from "vue";
import {jsPlumb} from "jsplumb";


const props = defineProps({graph: {type: Object}})
console.log(`[props]: `, )
const graph = props.graph
const tables = ref(graph.data.tables);
const relations = ref(graph.data.relations);
const state = {
  tables: graph.data.tables,
  relations: graph.data.relations,
};
const activeId = ref(false)
let instance;

onMounted(() => {
  nextTick(() => {
    console.log(`[jsPlumb]: `, jsPlumb);
    jsPlumb.ready(function () {
      jsPlumb.batch(() => {
        // 配置
        const connectorPaintStyle = {
          strokeWidth: 1,
          stroke: "#61B7CF",
          joinstyle: "round",
          outlineStroke: "white",
          outlineWidth: 0,
        };
        const connectorHoverStyle = {
          strokeWidth: 2,
          stroke: "#216477",
          outlineWidth: 0,
          outlineStroke: "white",
        };
        const endpointHoverStyle = {
          fill: "#216477",
          stroke: "#216477",
        };
        const sourceEndpoint = {
          endpoint: "Blank",
          paintStyle: {
            stroke: "#7AB02C",
            fill: "transparent",
            radius: 1,
            strokeWidth: 1,
          },
          isSource: true,
          connector: ["Bezier", {curviness: 180}],
          // connectorStyle: connectorPaintStyle,
          // hoverPaintStyle: endpointHoverStyle,
          // connectorHoverStyle: connectorHoverStyle,
          // maxConnections: -1,
          dragOptions: {},
          // overlays: [
          //   ["Label", {
          //     location: [0.5, 1.5],
          //     label: "",
          //     cssClass: "endpointSourceLabel",
          //     visible: true
          //   }]
          // ]
        };
        const targetEndpoint = {
          endpoint: "Blank",
          paintStyle: {fill: "#7AB02C", radius: 2},
          hoverPaintStyle: endpointHoverStyle,
          maxConnections: -1,
          dropOptions: {hoverClass: "hover", activeClass: "active"},
          isTarget: true,
          overlays: [
            [
              "Label",
              {
                location: [0.5, -0.5],
                label: "",
                cssClass: "endpointTargetLabel",
                visible: true,
              },
            ],
          ],
        };

        instance = jsPlumb.getInstance({
          DragOptions: {cursor: 'pointer', zIndex: 2000},
          ConnectionsDetachable: false,
          LogEnabled: true,
          Container: "container",
          Connector: ["Bezier", {curviness: 150}],
          Anchors: ["RightMiddle", "LeftMiddle"],
          // 绘制线
          PaintStyle: {
            stroke: "#ccc",
            strokeWidth: 0.1, //线的宽度
          },
          Endpoint: [
            "Blank",
            {
              height: 0,
              width: 0,
            },
          ],
          ConnectionOverlays: [
            //连线的叠加组件，如箭头、标签
            [
              "PlainArrow",
              {
                location: 1,
                visible: true,
                width: 9,
                length: 12,
                id: "ARROW",
                events: {
                  click: function () {
                    // alert("you clicked on the arrow overlay");
                  },
                },
              },
            ],
          ],
        });

        const fix = {
          paintStyle: {stroke: '#34495e', strokeWidth: 2},
          hoverPaintStyle: {stroke: '#34495e', strokeWidth: 2},
          connectorStyle: connectorHoverStyle,
          connectorHoverStyle: connectorHoverStyle
        };

        const hover = {
          paintStyle: {stroke: '#999', strokeWidth: 2},
          hoverPaintStyle: {stroke: '#999', strokeWidth: 2},
          connectorStyle: {
            strokeWidth: 2,
            stroke: "#999",
            outlineWidth: 0,
            outlineStroke: "white",
          },
          connectorHoverStyle: connectorHoverStyle
        };
        instance.registerConnectionType('fixed', fix);
        instance.registerConnectionType('hover', hover);
        instance.lineObject = {}
        instance.lineObject1 = {}
        window.instance = window.jsp = instance
        initTable(state.tables);


      })
    });
  });
});

// 配置表
function initTable(tables, sourceEndpoint = {}, targetEndpoint = {}) {
  for (const tablesKey in tables) {
    const table = tables[tablesKey];
    const {id} = table;
    initEndPoint(id, state.relations, sourceEndpoint, targetEndpoint);
  }
  instance.draggable(jsPlumb.getSelector(".table"), {grid: [20, 20]});
}

function initEndPoint(id, cols, sourceEndpoint, targetEndpoint) {
  for (let c of cols) {
    const {sourceColumn, sourceTable, targetColumn, targetTable} = c;
    const line = instance.connect({
      source: `${sourceTable}_${sourceColumn}`,
      target: `${targetTable}_${targetColumn}`,
    });

    instance.lineObject[`${sourceTable}_${sourceColumn}@${targetTable}_${targetColumn}`] = line
    instance.lineObject1[`${sourceTable}___${sourceColumn}@${targetTable}___${targetColumn}`] = line

    instance.addEndpoint(`${targetTable}_${targetColumn}`, targetEndpoint)
    instance.addEndpoint(`${sourceTable}_${sourceColumn}`, sourceEndpoint)
  }
}


const currentNode = reactive({
  target: [],
  source: []
})

const fixedSelectNode = ref([])

function fixedSelectHighlight() {
  currentNode.target.map(item => {
    const tid = `${item.targetTable}_${item.targetColumn}`
    const sid = `${item.sourceTable}_${item.sourceColumn}`
    highlight(sid + '@' + tid)
  })
  currentNode.source.map(item => {
    const tid = `${item.targetTable}_${item.targetColumn}`
    const sid = `${item.sourceTable}_${item.sourceColumn}`
    highlight(sid + '@' + tid)
  })
}


//Highlight
function highlight(lineObjectKey, type = 'fixed') {
  instance.lineObject[lineObjectKey]?.setType(type)
}  //Highlight
function dehighlight(lineObjectKey) {
  for (const lineObjectKey in instance.lineObject) {
    instance.lineObject[lineObjectKey].setType('default')
  }
}


function handleOver(e) {
  instance.batch(() => {
    dehighlight()
    const elId = e.target.id;
    const target = state.relations.map(
        (item) => `${item.sourceTable}_${item.sourceColumn}` === e.target.id
    );
    const source = state.relations.filter(
        (item) => `${item.sourceTable}_${item.sourceColumn}` === elId
    );

    const targets = findByDirection(elId, "target");
    const sources = findByDirection(elId, "source");

    // hover时也找全部的节点关系, 使用filter时只查找当前节点的上下一级
    targets
        .filter(item => {
          return item.sourceTable + '_' + item.sourceColumn === elId
        })
        .map(item => {
          const key = `${item.sourceTable}_${item.sourceColumn}@${item.targetTable}_${item.targetColumn}`
          highlight(key, 'hover')
        })

    sources
        .filter(item => {
          return item.targetTable + '_' + item.targetColumn === elId
        })
        .map(item => {
          const key = `${item.sourceTable}_${item.sourceColumn}@${item.targetTable}_${item.targetColumn}`
          highlight(key, 'hover')
        })

    fixedSelectHighlight()
  })
}

/**
 * 点击列名进行连线
 * @param e
 */
function handleClick(e) {
  dehighlight()
  const elId = e.target.id;
  activeId.value = elId
  document.querySelectorAll('.column-item').forEach(item => {
    item.classList.remove('active')
  })
  e.target.classList.add('active')
  const target = state.relations.map(
      (item) => `${item.sourceTable}_${item.sourceColumn}` === e.target.id
  );
  const source = state.relations.filter(
      (item) => `${item.sourceTable}_${item.sourceColumn}` === elId
  );

  const targets = findByDirection(elId, "target");
  const sources = findByDirection(elId, "source");


  targets.map(item => {
    const tid = `${item.targetTable}_${item.targetColumn}`
    const sid = `${item.sourceTable}_${item.sourceColumn}`
    highlight(sid + '@' + tid)
  })
  sources.map(item => {
    const tid = `${item.targetTable}_${item.targetColumn}`
    const sid = `${item.sourceTable}_${item.sourceColumn}`
    highlight(sid + '@' + tid)
  })
  currentNode.source = sources
  currentNode.target = targets

}

/**
 * 按节点方向查询
 * @param id
 * @param direction {'target' | 'source'}
 * @returns {*[]}
 */
function findByDirection(id, direction = "target") {
  const result = [];
  const dataMap = {};

  function _find(id) {
    const relations = state.relations.filter((item) => {
      dataMap.target = item.sourceTable + "_" + item.sourceColumn;
      dataMap.source = item.targetTable + "_" + item.targetColumn;
      if (dataMap[direction] === id) {
        return item;
      }
    });
    if (relations.length) {
      result.push(...relations);
      relations.map((item) =>
          direction === "target"
              ? _find(item.targetTable + "_" + item.targetColumn)
              : _find(item.sourceTable + "_" + item.sourceColumn)
      );
    }
  }

  _find(id);
  return result;
}

function handleSelect(event) {
  instance.batch(() => {
    fixedSelectNode.value.push(currentNode)
    fixedSelectHighlight()
  })
}
</script>

<template>
  <div id="container" class="container">
    <div
        v-for="(table, index) of state.tables"
        :id="table.id"
        :key="index"
        :style="{ left: table.left + 'px', top: table.top + 'px' }"
        class="table"
    >
      <div class="table-head">{{ table.id }}</div>

      <ul class="column-box">
        <li
            v-for="(col, i) of Object.keys(table.columns)"
            :id="`${table.id}_${col}`"
            :key="i"
            class="column-item"
            @click="handleClick"
            @mouseenter="handleOver"
            @drog.prevent
        >
          {{ col }}
        </li>
      </ul>
    </div>
  </div>
</template>

<style lang="scss" scoped>
.container {
  zoom: 75%;
  min-height: 750px;
  height: 750px;
  max-height: 1600px;
  border: 1px solid #CCC;
  background-color: var(--bg-color);
  /*background-color:var(--bg-grey);*/
  display: flex;
  position: relative;

  .table {
    width: 160px;
    position: absolute;
    display: inline-block;
    min-width: 100px;
    background-color: #27ae60;
    border: 1px solid #27ae60;
    -webkit-border-radius: 4px;
    -moz-border-radius: 4px;
    border-radius: 4px;
    transform: translate3d(0, 0, 0);
    backface-visibility: hidden;


    &-head {
      color: #fff;
      font-size: 18px;
      text-align: center;
    }

    .column-box {
      background-color: #fff;
      font-size: 14px;

      .column-item {
        padding: 4px;

        &:hover {
          background-color: skyblue;
        }
      }

    }
  }
}

.active {
  background-color: green !important;
  color: #ffffff;
}
</style>
