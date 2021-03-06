// r124

const assets = [
	'./',

	'./manifest.json',
	'./images/icon.png',

	'../files/favicon.ico',

	'../build/three.module.js',

	'../examples/jsm/controls/TransformControls.js',

	'../examples/jsm/libs/chevrotain.module.min.js',
	'../examples/jsm/libs/inflate.module.min.js',
	'../examples/jsm/libs/jszip.module.min.js',

	'../examples/js/libs/draco/draco_decoder.js',
	'../examples/js/libs/draco/draco_decoder.wasm',
	'../examples/js/libs/draco/draco_encoder.js',
	'../examples/js/libs/draco/draco_wasm_wrapper.js',

	'../examples/js/libs/draco/gltf/draco_decoder.js',
	'../examples/js/libs/draco/gltf/draco_decoder.wasm',
	'../examples/js/libs/draco/gltf/draco_wasm_wrapper.js',

	'../examples/jsm/libs/rhino3dm/rhino3dm.wasm',
	'../examples/jsm/libs/rhino3dm/rhino3dm.js',

	'../examples/jsm/loaders/3DMLoader.js',
	'../examples/jsm/loaders/3MFLoader.js',
	'../examples/jsm/loaders/AMFLoader.js',
	'../examples/jsm/loaders/ColladaLoader.js',
	'../examples/jsm/loaders/DRACOLoader.js',
	'../examples/jsm/loaders/FBXLoader.js',
	'../examples/jsm/loaders/GLTFLoader.js',
	'../examples/jsm/loaders/KMZLoader.js',
	'../examples/jsm/loaders/MD2Loader.js',
	'../examples/jsm/loaders/OBJLoader.js',
	'../examples/jsm/loaders/MTLLoader.js',
	'../examples/jsm/loaders/PLYLoader.js',
	'../examples/jsm/loaders/RGBELoader.js',
	'../examples/jsm/loaders/STLLoader.js',
	'../examples/jsm/loaders/SVGLoader.js',
	'../examples/jsm/loaders/TGALoader.js',
	'../examples/jsm/loaders/TDSLoader.js',
	'../examples/jsm/loaders/VOXLoader.js',
	'../examples/jsm/loaders/VRMLLoader.js',
	'../examples/jsm/loaders/VTKLoader.js',
	'../examples/jsm/loaders/XYZLoader.js',

	'../examples/jsm/curves/NURBSCurve.js',
	'../examples/jsm/curves/NURBSUtils.js',

	'../examples/jsm/exporters/ColladaExporter.js',
	'../examples/jsm/exporters/DRACOExporter.js',
	'../examples/jsm/exporters/GLTFExporter.js',
	'../examples/jsm/exporters/OBJExporter.js',
	'../examples/jsm/exporters/PLYExporter.js',
	'../examples/jsm/exporters/STLExporter.js',

	'../examples/jsm/helpers/VertexNormalsHelper.js',

	'../examples/jsm/geometries/TeapotBufferGeometry.js',

	'../examples/jsm/webxr/VRButton.js',

	'./images/rotate.svg',
	'./images/scale.svg',
	'./images/translate.svg',

	'./js/libs/codemirror/codemirror.css',
	'./js/libs/codemirror/theme/monokai.css',

	'./js/libs/codemirror/codemirror.js',
	'./js/libs/codemirror/mode/javascript.js',
	'./js/libs/codemirror/mode/glsl.js',

	'./js/libs/esprima.js',
	'./js/libs/jsonlint.js',

	'./js/libs/codemirror/addon/dialog.css',
	'./js/libs/codemirror/addon/show-hint.css',
	'./js/libs/codemirror/addon/tern.css',

	'./js/libs/codemirror/addon/dialog.js',
	'./js/libs/codemirror/addon/show-hint.js',
	'./js/libs/codemirror/addon/tern.js',
	'./js/libs/acorn/acorn.js',
	'./js/libs/acorn/acorn_loose.js',
	'./js/libs/acorn/walk.js',
	'./js/libs/ternjs/polyfill.js',
	'./js/libs/ternjs/signal.js',
	'./js/libs/ternjs/tern.js',
	'./js/libs/ternjs/def.js',
	'./js/libs/ternjs/comment.js',
	'./js/libs/ternjs/infer.js',
	'./js/libs/ternjs/doc_comment.js',
	'./js/libs/tern-threejs/threejs.js',

	'./js/libs/signals.min.js',
	'./js/libs/ui.js',
	'./js/libs/ui.three.js',
	'./js/libs/reconnecting-websocket.js',

	'./js/libs/app.js',
	'./js/Player.js',
	'./js/Script.js',

	//

	'./css/main.css',

	'./js/EditorControls.js',
	'./js/Storage.js',

	'./js/Editor.js',
	'./js/Config.js',
	'./js/History.js',
	'./js/Loader.js',
	'./js/LoaderUtils.js',
	'./js/BasiliskBufferGeometry.js',
	'./js/Menubar.js',
	'./js/Menubar.Basilisk.js',
	'./js/Menubar.File.js',
	'./js/Menubar.Edit.js',
	'./js/Menubar.Add.js',
	'./js/Menubar.Play.js',
	'./js/Menubar.Examples.js',
	'./js/Menubar.Help.js',
	'./js/Menubar.Status.js',
	'./js/Resizer.js',
	'./js/Sidebar.js',
	'./js/Sidebar.Basilisk.js',
	'./js/Sidebar.Scene.js',
	'./js/Sidebar.Project.js',
	'./js/Sidebar.Settings.js',
	'./js/Sidebar.Settings.Shortcuts.js',
	'./js/Sidebar.Settings.Viewport.js',
	'./js/Sidebar.Properties.js',
	'./js/Sidebar.Object.js',
	'./js/Sidebar.Geometry.js',
	'./js/Sidebar.Geometry.Geometry.js',
	'./js/Sidebar.Geometry.BufferGeometry.js',
	'./js/Sidebar.Geometry.BasiliskBufferGeometry.js',
	'./js/Sidebar.Geometry.Modifiers.js',
	'./js/Sidebar.Geometry.BoxBufferGeometry.js',
	'./js/Sidebar.Geometry.CircleBufferGeometry.js',
	'./js/Sidebar.Geometry.CylinderBufferGeometry.js',
	'./js/Sidebar.Geometry.DodecahedronBufferGeometry.js',
	'./js/Sidebar.Geometry.ExtrudeBufferGeometry.js',
	'./js/Sidebar.Geometry.IcosahedronBufferGeometry.js',
	'./js/Sidebar.Geometry.LatheBufferGeometry.js',
	'./js/Sidebar.Geometry.OctahedronBufferGeometry.js',
	'./js/Sidebar.Geometry.PlaneBufferGeometry.js',
	'./js/Sidebar.Geometry.RingBufferGeometry.js',
	'./js/Sidebar.Geometry.SphereBufferGeometry.js',
	'./js/Sidebar.Geometry.ShapeBufferGeometry.js',
	'./js/Sidebar.Geometry.TetrahedronBufferGeometry.js',
	'./js/Sidebar.Geometry.TorusBufferGeometry.js',
	'./js/Sidebar.Geometry.TorusKnotBufferGeometry.js',
	'./js/Sidebar.Geometry.TubeBufferGeometry.js',
	'./js/Sidebar.Geometry.TeapotBufferGeometry.js',
	'./js/Sidebar.Material.js',
	'./js/Sidebar.Animation.js',
	'./js/Sidebar.Script.js',
	'./js/Sidebar.History.js',
	'./js/Strings.js',
	'./js/Toolbar.js',
	'./js/Viewport.js',
	'./js/Viewport.Camera.js',
	'./js/Viewport.Info.js',
	'./js/Viewport.ViewHelper.js',

	'./js/Command.js',
	'./js/commands/AddObjectCommand.js',
	'./js/commands/RemoveObjectCommand.js',
	'./js/commands/MoveObjectCommand.js',
	'./js/commands/SetPositionCommand.js',
	'./js/commands/SetRotationCommand.js',
	'./js/commands/SetScaleCommand.js',
	'./js/commands/SetValueCommand.js',
	'./js/commands/SetUuidCommand.js',
	'./js/commands/SetColorCommand.js',
	'./js/commands/SetGeometryCommand.js',
	'./js/commands/SetGeometryValueCommand.js',
	'./js/commands/MultiCmdsCommand.js',
	'./js/commands/AddScriptCommand.js',
	'./js/commands/RemoveScriptCommand.js',
	'./js/commands/SetScriptValueCommand.js',
	'./js/commands/SetMaterialCommand.js',
	'./js/commands/SetMaterialColorCommand.js',
	'./js/commands/SetMaterialMapCommand.js',
	'./js/commands/SetMaterialValueCommand.js',
	'./js/commands/SetMaterialVectorCommand.js',
	'./js/commands/SetSceneCommand.js',
	'./js/commands/Commands.js',

	//

	'./examples/arkanoid.app.json',
	'./examples/camera.app.json',
	'./examples/particles.app.json',
	'./examples/pong.app.json',
	'./examples/shaders.app.json'

];

self.addEventListener( 'install', async function () {

	const cache = await caches.open( 'threejs-editor' );

	assets.forEach( function ( asset ) {

		cache.add( asset ).catch( function () {

			console.error( '[SW] Cound\'t cache:', asset );

		} );

	} );

} );

self.addEventListener( 'fetch', async function ( event ) {

	const request = event.request;
	event.respondWith( cacheFirst( request ) );

} );

async function cacheFirst( request ) {

	const cachedResponse = await caches.match( request );

	if ( cachedResponse === undefined ) {

		console.error( '[SW] Not cached:', request.url );
		return fetch( request );

	}

	return cachedResponse;

}
