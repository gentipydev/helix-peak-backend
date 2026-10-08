// A protein resolved on demand, read from the live service through the app's
// production client, then walked on the real screen. It is kept here, in the
// backend's scripts/, and never committed to the app. From helix-peek/:
//
//   cp ../helix-peek-backend/scripts/live_resolved_check.dart test/
//   LIVE_BACKEND=https://helix-peak-backend.onrender.com LIVE_RESOLVED=B2M \
//     flutter test test/live_resolved_check.dart
//   rm test/live_resolved_check.dart
//
// LIVE_EXPECT=refused is for a protein whose ESM-2 track the scorer declined.
import 'dart:io';

import 'package:flutter/material.dart';
import 'package:flutter_test/flutter_test.dart';
import 'package:helixpeek/core/biology/gene_record.dart';
import 'package:helixpeek/core/catalog/protein_resolver.dart';
import 'package:helixpeek/core/catalog/protein_suggestion.dart';
import 'package:helixpeek/core/catalog/protein_target.dart';
import 'package:helixpeek/core/catalog/protein_track.dart';
import 'package:helixpeek/core/evidence/protein_constraint.dart';
import 'package:helixpeek/core/network/dio_api_client.dart';
import 'package:helixpeek/core/network/track_client.dart';
import 'package:helixpeek/core/theme/app_theme.dart';
import 'package:helixpeek/features/gene_lookup/data/datasources/gene_remote_data_source.dart';
import 'package:helixpeek/features/gene_lookup/data/models/gene_record_dto.dart';
import 'package:helixpeek/features/gene_lookup/data/repositories/protein_catalog_repository.dart';
import 'package:helixpeek/features/gene_lookup/presentation/anatomy/anatomy_canvas.dart';
import 'package:helixpeek/features/gene_lookup/presentation/anatomy/anatomy_screen.dart';
import 'package:helixpeek/features/gene_lookup/presentation/constraint/constraint_toolbar.dart';
import 'package:helixpeek/shared/anatomy/anatomy_stages.dart';

const Size _phone = Size(390, 844);

void main() {
  final String baseUrl = Platform.environment['LIVE_BACKEND'] ?? '';
  final String gene = Platform.environment['LIVE_RESOLVED'] ?? '';
  final bool scored = (Platform.environment['LIVE_EXPECT'] ?? 'ready') == 'ready';
  final bool on = baseUrl.isNotEmpty && gene.isNotEmpty;
  final String slug = gene.toLowerCase();
  final TrackState expected = scored ? TrackState.ready : TrackState.refused;

  late Directory root;
  late ProteinTarget target;
  late GeneRecord record;
  ProteinConstraint? constraint;
  Object? constraintError;
  late ResolveStatus status;
  late TrackState trackState;
  late SuggestionPage suggested;
  late int listed;
  late bool inList;

  setUpAll(() async {
    if (!on) return;
    // The widget binding answers every request with a 400; this is the one
    // file that wants the real network.
    HttpOverrides.global = null;
    root = Directory.systemTemp.createTempSync('helixpeek-resolved');
    final DioApiClient api = DioApiClient(
      baseUrl: baseUrl,
      timeout: const Duration(seconds: 60),
    );
    final ProteinCatalogRepository catalog = ProteinCatalogRepository(api, null);
    await catalog.refresh();
    listed = catalog.all.length;
    inList = catalog.bySlug(slug) != null;

    // The row, as a deep link outside the catalog reads it.
    target = await catalog.protein(slug);

    // The record and the ESM-2 track, as the walk reads them.
    final TrackClient tracks = TrackClient(
      api,
      cache: Directory('${root.path}/tracks')..createSync(),
    );
    record = (await TrackGeneDataSource(tracks).fetchGene(target.query))
        .toEntity();
    try {
      constraint = await ProteinConstraint.load(target, tracks: tracks);
    } on Object catch (error) {
      constraintError = error;
    }

    // What the search screen asks.
    final ProteinResolver resolver = ProteinResolver(api);
    status = await resolver.status(gene);
    trackState = await resolver.trackState(slug, TrackKind.constraint);
    suggested = await resolver.suggest(gene);
    catalog.dispose();
  });

  tearDownAll(() {
    if (on) root.deleteSync(recursive: true);
  });

  test('the service serves it as a built protein outside the list', () {
    if (!on) {
      markTestSkipped('set LIVE_BACKEND and LIVE_RESOLVED');
      return;
    }
    expect(listed, 20);
    expect(inList, isFalse);
    expect((target.slug, target.gene), (slug, gene));
    expect(target.summary, contains('Built on demand from UniProt'));
    expect(target.state(TrackKind.record), TrackState.ready);
    expect(target.state(TrackKind.constraint), expected);
    expect(target.state(TrackKind.structure), TrackState.absent);
    expect(target.state(TrackKind.clinvar), TrackState.absent);
    expect(target.scored, scored);

    expect(record.gene, gene);
    expect(record.protein!.translation.length, target.facts.residues);
    if (scored) {
      expect(constraintError, isNull);
      expect(constraint!.positions.length, target.facts.residues);
    } else {
      // ignore: avoid_print
      print('$gene: asking for its ESM-2 track gives: $constraintError');
      expect(constraint, isNull);
    }

    expect(status.state, ResolveState.ready);
    expect(status.slug, slug);
    expect(trackState, expected);
    final ProteinSuggestion row = suggested.suggestions.firstWhere(
      (ProteinSuggestion s) => s.gene == gene,
    );
    expect(row.status, SuggestionStatus.ready);
    expect(row.slug, slug);
    // ignore: avoid_print
    print(
      '$gene: ${target.display}, ${target.facts.residues} residues, '
      '${record.exons.length} exons, ${target.facts.chains} chain(s), '
      '${target.facts.bridges} bridge(s); record ${record.lengthBp} bp'
      '${record.isIntronCompressed ? ' (introns shortened)' : ''}; '
      'ESM-2 track ${expected.name}.',
    );
  });

  testWidgets('the walk draws every page of it', (WidgetTester tester) async {
    if (!on) {
      markTestSkipped('set LIVE_BACKEND and LIVE_RESOLVED');
      return;
    }
    addTearDown(() => tester.binding.setSurfaceSize(null));
    await tester.binding.setSurfaceSize(_phone);
    await tester.pumpWidget(
      MaterialApp(
        theme: AppTheme.analysis,
        debugShowCheckedModeBanner: false,
        home: MediaQuery(
          data: const MediaQueryData(disableAnimations: true),
          child: AnatomyScreen(
            target: target,
            record: record,
            constraint: constraint,
          ),
        ),
      ),
    );
    await tester.pumpAndSettle();
    expect(find.byType(AnatomyCanvas), findsOneWidget);

    final List<AnatomyStage> stages = AnatomyModel.derive(
      record,
      chain: target.chain,
    ).stages;
    final int pages = stages.length + 1;
    final int proteinPage = stages.indexWhere(
      (AnatomyStage s) => s.kind == StageKind.protein,
    );
    expect(proteinPage, isNonNegative);
    final Set<String> seen = <String>{};
    for (final bool forward in <bool>[true, false]) {
      for (int i = 1; i < pages; i++) {
        final Rect screen = tester.getRect(find.byType(AnatomyScreen));
        await tester.dragFrom(
          Offset(screen.center.dx, screen.bottom - 40),
          Offset(forward ? -160 : 160, 0),
        );
        await tester.pumpAndSettle();
        expect(tester.takeException(), isNull, reason: 'page $i');
        if (forward && i == proteinPage) {
          expect(
            find.byType(ConstraintToolbar),
            scored ? findsOneWidget : findsNothing,
          );
          // What the protein page says in words, to read what a reader is
          // told about a track that is not there.
          for (final Element e in find.byType(Text).evaluate()) {
            final String? words = (e.widget as Text).data;
            if (words != null && RegExp('ESM|scor|conserv', caseSensitive: false).hasMatch(words)) {
              seen.add(words);
            }
          }
        }
      }
    }
    // ignore: avoid_print
    print('$gene: walked $pages pages forward and back. Protein page says: $seen');
  });
}
