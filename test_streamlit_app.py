import csv
import io
import unittest
from pathlib import Path

import numpy as np
from streamlit.testing.v1 import AppTest
from streamlit_app import calculate_model, preview_figure, results_csv
from TPMS_Analyzer import analyze


class WebAppTests(unittest.TestCase):
    def test_results_match_desktop_and_sphere_is_fixed(self):
        model = calculate_model('Solid', 'Gyroid', 70., '2.54', 24, 'auto')
        expected = analyze('Solid', 'Gyroid', 70., 2.54, 24, 4)
        self.assertEqual(model['rows'], expected.rows)
        self.assertNotIn('field', model)
        figures = [preview_figure(model, state) for state in (True, False)]
        spheres = [next(t for t in f.data if t.name == 'Representative pore') for f in figures]
        for coordinate in ('x', 'y', 'z'):
            np.testing.assert_array_equal(getattr(spheres[0], coordinate),
                                          getattr(spheres[1], coordinate))
        self.assertEqual(spheres[0].opacity, 1.)
        self.assertLess(figures[0].layout.scene.xaxis.range[0], 0)
        self.assertGreater(figures[0].layout.scene.xaxis.range[1], model['alpha'])
        glossy = preview_figure(model, False, 'Glossy')
        self.assertGreater(glossy.data[0].lighting.specular, figures[0].data[0].lighting.specular)
        rows = list(csv.reader(io.StringIO(results_csv(model['rows']))))
        self.assertEqual(len(rows), 18)

    def test_calculate_toggle_and_invalid_input(self):
        app = AppTest.from_file(str(Path(__file__).with_name('streamlit_app.py'))).run()
        self.assertFalse(app.exception)
        self.assertTrue(app.info)
        app.number_input[1].set_value(24)
        app.button[0].click().run(timeout=30)
        self.assertFalse(app.exception)
        self.assertEqual(len(app.metric), 3)
        before = app.session_state['analysis']
        app.toggle[0].set_value(False).run()
        self.assertFalse(app.exception)
        np.testing.assert_array_equal(before['vertices'], app.session_state['analysis']['vertices'])
        self.assertEqual(before['rows'], app.session_state['analysis']['rows'])
        next(s for s in app.selectbox if s.label == 'Surface finish').set_value('Glossy').run()
        self.assertFalse(app.exception)
        self.assertEqual(before['rows'], app.session_state['analysis']['rows'])
        app.text_input[0].set_value('invalid')
        app.button[0].click().run()
        self.assertTrue(app.error)
        self.assertFalse(app.exception)
        self.assertEqual(before['rows'], app.session_state['analysis']['rows'])


if __name__ == '__main__':
    unittest.main()
