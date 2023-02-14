<script setup>
import {nextTick, onMounted, reactive, ref} from "vue";
import {jsPlumb} from "jsplumb";

var graph = {
  data: {
    tables: {
      RESULT__T1: {
        id: "T15",
        alias: "",
        name: "RESULT__T1",
        columns: {
          A: "VARCHAR",
          B: "VARCHAR",
          CODE: "VARCHAR",
          C: "VARCHAR",
          D: "VARCHAR",
          ID: "VARCHAR",
        },
        type: "o",
        top: 100,
        left: 1000,
        level: 2,
        groups: [],
        targets: ["T18"],
        metaColumns: {},
        parent: "RESULT",
      },
      RESULT__T2: {
        id: "T16",
        alias: "",
        name: "RESULT__T2",
        columns: {
          CODE: "VARCHAR",
          C: "VARCHAR",
          ID: "VARCHAR",
        },
        type: "o",
        top: 320,
        left: 1000,
        level: 2,
        groups: [],
        targets: ["T18"],
        metaColumns: {},
        parent: "RESULT",
      },
      TAB_2: {
        id: "T11",
        alias: "",
        name: "TAB_2",
        columns: {
          CODE: "VARCHAR",
          C: "VARCHAR",
          ID: "VARCHAR",
        },
        type: "g",
        top: 520,
        left: 600,
        level: 1,
        groups: [],
        targets: ["T16"],
        metaColumns: {},
        parent: "RESULT__T2",
      },
      TAB_3: {
        id: "T12",
        alias: "",
        name: "TAB_3",
        columns: {
          ID: "VARCHAR",
          NAME: "VARCHAR",
        },
        type: "g",
        top: 480,
        left: 1000,
        level: 2,
        groups: [],
        targets: ["T18"],
        metaColumns: {},
        parent: "RESULT",
      },
      TAB_0: {
        id: "T13",
        alias: "",
        name: "TAB_0",
        columns: {
          AA: "VARCHAR",
          COM: "VARCHAR",
          B: "VARCHAR",
          C: "VARCHAR",
          ID: "VARCHAR",
          DESC: "VARCHAR",
        },
        type: "g",
        top: 190,
        left: 200,
        level: 0,
        groups: [],
        targets: ["T17"],
        metaColumns: {},
        parent: "RESULT__T1__U2",
      },
      TAB_1: {
        id: "T14",
        alias: "",
        name: "TAB_1",
        columns: {
          A: "VARCHAR",
          B: "VARCHAR",
          C: "VARCHAR",
          D: "VARCHAR",
          ID: "VARCHAR",
        },
        type: "g",
        top: 410,
        left: 200,
        level: 0,
        groups: [],
        targets: ["T19"],
        metaColumns: {},
        parent: "RESULT__T1__U1",
      },
      RESULT__T1__U2: {
        id: "T17",
        alias: "",
        name: "RESULT__T1__U2",
        columns: {
          A: "VARCHAR",
          B: "VARCHAR",
          CODE: "VARCHAR",
          C: "VARCHAR",
          D: "VARCHAR",
          ID: "VARCHAR",
        },
        type: "o",
        top: 80,
        left: 600,
        level: 1,
        groups: [],
        targets: ["T15"],
        metaColumns: {},
        parent: "RESULT__T1",
      },
      RESULT: {
        id: "T18",
        alias: "",
        name: "RESULT",
        columns: {
          AAA: "VARCHAR",
          CCC: "VARCHAR",
          BBB: "VARCHAR",
          EEE: "VARCHAR",
          DDD: "VARCHAR",
          FFF: "VARCHAR",
        },
        type: "b",
        top: 240,
        left: 1500,
        level: 3,
        groups: [],
        targets: [],
        metaColumns: {},
        parent: null,
      },
      RESULT__T1__U1: {
        id: "T19",
        alias: "",
        name: "RESULT__T1__U1",
        columns: {
          A: "VARCHAR",
          B: "VARCHAR",
          CODE: "VARCHAR",
          C: "VARCHAR",
          D: "VARCHAR",
          ID: "VARCHAR",
        },
        type: "o",
        top: 300,
        left: 600,
        level: 1,
        groups: [],
        targets: ["T15"],
        metaColumns: {},
        parent: "RESULT__T1",
      },
    },
    relations: [
      {
        sourceTable: "T19",
        targetTable: "T15",
        sourceColumn: "A",
        targetColumn: "A",
        type: "A",
      },
      {
        sourceTable: "T19",
        targetTable: "T15",
        sourceColumn: "B",
        targetColumn: "B",
        type: "A",
      },
      {
        sourceTable: "T19",
        targetTable: "T15",
        sourceColumn: "C",
        targetColumn: "C",
        type: "A",
      },
      {
        sourceTable: "T19",
        targetTable: "T15",
        sourceColumn: "D",
        targetColumn: "D",
        type: "A",
      },
      {
        sourceTable: "T19",
        targetTable: "T15",
        sourceColumn: "ID",
        targetColumn: "ID",
        type: "A",
      },
      {
        sourceTable: "T19",
        targetTable: "T15",
        sourceColumn: "CODE",
        targetColumn: "CODE",
        type: "A",
      },
      {
        sourceTable: "T17",
        targetTable: "T15",
        sourceColumn: "A",
        targetColumn: "A",
        type: "A",
      },
      {
        sourceTable: "T17",
        targetTable: "T15",
        sourceColumn: "B",
        targetColumn: "B",
        type: "A",
      },
      {
        sourceTable: "T17",
        targetTable: "T15",
        sourceColumn: "C",
        targetColumn: "C",
        type: "A",
      },
      {
        sourceTable: "T17",
        targetTable: "T15",
        sourceColumn: "D",
        targetColumn: "D",
        type: "A",
      },
      {
        sourceTable: "T17",
        targetTable: "T15",
        sourceColumn: "ID",
        targetColumn: "ID",
        type: "A",
      },
      {
        sourceTable: "T17",
        targetTable: "T15",
        sourceColumn: "CODE",
        targetColumn: "CODE",
        type: "A",
      },
      {
        sourceTable: "T11",
        targetTable: "T16",
        sourceColumn: "C",
        targetColumn: "C",
        type: "A",
      },
      {
        sourceTable: "T11",
        targetTable: "T16",
        sourceColumn: "CODE",
        targetColumn: "CODE",
        type: "A",
      },
      {
        sourceTable: "T11",
        targetTable: "T16",
        sourceColumn: "ID",
        targetColumn: "ID",
        type: "A",
      },
      {
        sourceTable: "T13",
        targetTable: "T17",
        sourceColumn: "AA",
        targetColumn: "A",
        type: "A",
      },
      {
        sourceTable: "T13",
        targetTable: "T17",
        sourceColumn: "B",
        targetColumn: "B",
        type: "A",
      },
      {
        sourceTable: "T13",
        targetTable: "T17",
        sourceColumn: "C",
        targetColumn: "C",
        type: "A",
      },
      {
        sourceTable: "T13",
        targetTable: "T17",
        sourceColumn: "COM",
        targetColumn: "D",
        type: "A",
      },
      {
        sourceTable: "T13",
        targetTable: "T17",
        sourceColumn: "ID",
        targetColumn: "ID",
        type: "A",
      },
      {
        sourceTable: "T13",
        targetTable: "T17",
        sourceColumn: "DESC",
        targetColumn: "CODE",
        type: "A",
      },
      {
        sourceTable: "T15",
        targetTable: "T18",
        sourceColumn: "A",
        targetColumn: "AAA",
        type: "A",
      },
      {
        sourceTable: "T15",
        targetTable: "T18",
        sourceColumn: "A",
        targetColumn: "BBB",
        type: "A",
      },
      {
        sourceTable: "T15",
        targetTable: "T18",
        sourceColumn: "B",
        targetColumn: "BBB",
        type: "A",
      },
      {
        sourceTable: "T16",
        targetTable: "T18",
        sourceColumn: "C",
        targetColumn: "CCC",
        type: "A",
      },
      {
        sourceTable: "T15",
        targetTable: "T18",
        sourceColumn: "A",
        targetColumn: "DDD",
        type: "A",
      },
      {
        sourceTable: "T15",
        targetTable: "T18",
        sourceColumn: "A",
        targetColumn: "DDD",
        type: "A",
      },
      {
        sourceTable: "T15",
        targetTable: "T18",
        sourceColumn: "C",
        targetColumn: "EEE",
        type: "A",
      },
      {
        sourceTable: "T15",
        targetTable: "T18",
        sourceColumn: "A",
        targetColumn: "FFF",
        type: "A",
      },
      {
        sourceTable: "T12",
        targetTable: "T18",
        sourceColumn: "NAME",
        targetColumn: "FFF",
        type: "A",
      },
      {
        sourceTable: "T16",
        targetTable: "T18",
        sourceColumn: "CODE",
        targetColumn: "FFF",
        type: "A",
      },
      {
        sourceTable: "T14",
        targetTable: "T19",
        sourceColumn: "A",
        targetColumn: "A",
        type: "A",
      },
      {
        sourceTable: "T14",
        targetTable: "T19",
        sourceColumn: "B",
        targetColumn: "B",
        type: "A",
      },
      {
        sourceTable: "T14",
        targetTable: "T19",
        sourceColumn: "C",
        targetColumn: "C",
        type: "A",
      },
      {
        sourceTable: "T14",
        targetTable: "T19",
        sourceColumn: "D",
        targetColumn: "D",
        type: "A",
      },
      {
        sourceTable: "T14",
        targetTable: "T19",
        sourceColumn: "ID",
        targetColumn: "ID",
        type: "A",
      },
      {
        sourceTable: "T14",
        targetTable: "T19",
        sourceColumn: "ID",
        targetColumn: "CODE",
        type: "A",
      },
    ],
    tabRelations: [
      {
        sourceTable: "T15",
        targetTable: "T12",
        joinType: null,
        joinColumns: {
          ID: "ID",
        },
        joinColumnsStr: "{ID=ID}",
      },
      {
        sourceTable: "T15",
        targetTable: "T16",
        joinType: null,
        joinColumns: {
          CODE: "CODE",
          ID: "ID",
        },
        joinColumnsStr: "{CODE=CODE, ID=ID}",
      },
    ],
    TABID_DICT: {
      T12: "TAB_3",
      T11: "TAB_2",
      T14: "TAB_1",
      T13: "TAB_0",
      T16: "RESULT__T2",
      T15: "RESULT__T1",
      T18: "RESULT",
      T17: "RESULT__T1__U2",
      T19: "RESULT__T1__U1",
    },
  },
  message: null,
  code: 0,
};

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

        instance = window.jsp = jsPlumb.getInstance({
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
        console.log(`[jsp]: `, jsp)

        var fix = {
          paintStyle: {stroke: '#34495e', strokeWidth: 2},
          hoverPaintStyle: {stroke: '#34495e', strokeWidth: 2},
          connectorStyle: connectorHoverStyle,
          connectorHoverStyle: connectorHoverStyle
        };

        var hover = {
          paintStyle: {stroke: '#f00', strokeWidth: 2},
          hoverPaintStyle: {stroke: '#f00', strokeWidth: 2},
          connectorStyle: {
            strokeWidth: 2,
            stroke: "#f00",
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
        jsPlumb.fire("jsPlumbDemoLoaded", instance);
        instance.draggable(jsPlumb.getSelector(".table"), { grid: [20, 20]});



        function initEndPoint(id, cols) {
          for (let c of cols) {
            const {sourceColumn, sourceTable, targetColumn, targetTable} = c;

            instance.addEndpoint(`${id}_${c.targetColumn}`, {}, {
              anchor: 'RightMiddle', uuid: `${id}_${c.targetColumn}_RightMiddle`
            })
            instance.addEndpoint(`${id}_${c.sourceColumn}`, {}, {
              anchor: 'LeftMiddle', uuid: `${id}_${c.sourceColumn}_LeftMiddle`
            })


            const line = instance.connect({
              uuids: [
                  `${sourceTable}_${sourceColumn}_RightMiddle`,
                  `${targetTable}_${targetColumn}_LeftMiddle`
              ]
            })
            instance.lineObject[`${sourceTable}_${sourceColumn}@${targetTable}_${targetColumn}`] = line
            instance.lineObject1[`${sourceTable}___${sourceColumn}@${targetTable}___${targetColumn}`] = line

            // instance.addEndpoint(`${targetTable}_${targetColumn}`, targetEndpoint)
            // instance.addEndpoint(`${sourceTable}_${sourceColumn}`, sourceEndpoint)

            // instance.addEndpoint(`${targetTable}_${targetColumn}`, sourceEndpoint, {
            //   anchor: "RightMiddle", uuid: id + '_' + cols[c] + "_RightMiddle"
            // });
            // instance.addEndpoint(`${sourceTable}_${sourceColumn}}`, targetEndpoint, {
            //   anchor: "LeftMiddle", uuid: `${sourceTable}_${sourceColumn}` + "_LeftMiddle"
            // });
          }
        }
      })
    });
  });
});

// 配置表
function initTable(tables, sourceEndpoint = {}, targetEndpoint = {}) {

  for (const tablesKey in tables) {
    const table = tables[tablesKey];
    const {id} = table;

    // initEndPoint(id, state.relations);
    initAnchor(id, table);
    console.log(`[table]: `, table)
    initRelation(id, state.relations)
  }
}

function initAnchor(tableId, table) {
  for (const columnsKey in table.columns) {
    const elSelector = `${tableId}_${columnsKey}`
    instance.addEndpoint(elSelector, {}, {
      anchor: "RightMiddle",
      uuid: `${elSelector}_RightMiddle`
    })
    instance.addEndpoint(elSelector, {}, {
      anchor: "LeftMiddle",
      uuid: `${elSelector}_LeftMiddle`
    })
  }
}

function initRelation(tableId, relations) {
  relations.map(item => {
    console.log(`[item]: `, item)
    const {sourceTable, sourceColumn, targetTable, targetColumn} = item
    instance.connect({
      source: `${sourceTable}_${sourceColumn}`,
      target: `${targetTable}_${targetColumn}`,
      // anchors:["Right", "Left" ],
      // endpoint:"Rectangle",
    })
  })
}
function initPoint() {
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
        // .filter(item => {
        //   return item.sourceTable + '_' + item.sourceColumn === elId
        // })
        .map(item => {
          const key = `${item.sourceTable}_${item.sourceColumn}@${item.targetTable}_${item.targetColumn}`
          highlight(key, 'hover')
        })

    sources
        // .filter(item => {
        //   return item.targetTable + '_' + item.targetColumn === elId
        // })
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
            :class="['column-item', {active: `${table.id}_${col}` === false}]"
            @click="handleClick"
            @mouseover="handleOver"
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
