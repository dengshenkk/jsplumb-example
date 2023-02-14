import {jsPlumb} from "jsplumb";
import {nextTick, onMounted, reactive, ref} from "vue";


export function useJsplumb(container, data, relations) {

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


    const DEFAULT_SETTING = {
        DragOptions: {cursor: 'pointer', zIndex: 2000},
        ConnectionsDetachable: false,
        LogEnabled: true,
        Container: "",
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
    }

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
    }

    onMounted(() => {
        nextTick().then(_ => {
            DEFAULT_SETTING.Container = container
            const instance = jsPlumb.getInstance(DEFAULT_SETTING)
            // config
            instance.registerConnectionType('fixed', fix);
            instance.registerConnectionType('hover', hover);
            instance.lineObject = {}
            instance.lineObject1 = {}
            instance.ready(() => {
                // render
                instance.batch(() => {
                    initTable()
                })
            })
        })
    })
}
